#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

// [[Rcpp::plugins(cpp11)]]

namespace {

double projective_angle(
    const double perpendicular_x,
    const double perpendicular_y,
    const double pi) {
  double angle = std::atan2(perpendicular_y, perpendicular_x);
  angle = std::fmod(angle, pi);

  if (angle < 0.0) {
    angle += pi;
  }

  if (angle >= pi) {
    angle -= pi;
  }

  return angle;
}

} // namespace


// Construct the complete angular arrangement in dimension two.

// [[Rcpp::export]]
Rcpp::List angular_arrangement_d2_cpp(
   const Rcpp::NumericMatrix& Y,
   const double tol = 1e-12) {
   const int n = Y.nrow();
   const int d = Y.ncol();
   const double pi = std::acos(-1.0);

   if (d != 2) {
     Rcpp::stop("Y must have exactly two columns.");
   }

   if (n < 2) {
     Rcpp::stop("Y must contain at least two observations.");
   }

   if (!std::isfinite(tol) || tol < 0.0) {
     Rcpp::stop("tol must be a finite nonnegative number.");
   }

   double data_scale = 1.0;

   for (int row = 0; row < n; ++row) {
     for (int column = 0; column < 2; ++column) {
       const double value = Y(row, column);

       if (!std::isfinite(value)) {
         Rcpp::stop("Y must contain only finite values.");
       }

       data_scale = std::max(data_scale, std::abs(value));
     }
   }

   const double zero_threshold = tol * data_scale;

   std::vector<double> critical_angles;
   critical_angles.reserve(
     static_cast<std::size_t>(n) *
       static_cast<std::size_t>(n - 1)
   );

   // For every i < j, generate the two critical hyperplanes corresponding
   // to y_i + y_j and y_i - y_j.
   for (int i = 0; i < n - 1; ++i) {
     for (int j = i + 1; j < n; ++j) {
       const double sum_x = Y(i, 0) + Y(j, 0);
       const double sum_y = Y(i, 1) + Y(j, 1);

       const double difference_x = Y(i, 0) - Y(j, 0);
       const double difference_y = Y(i, 1) - Y(j, 1);

       const double sum_norm = std::hypot(sum_x, sum_y);

       if (sum_norm > zero_threshold) {
         critical_angles.push_back(
           projective_angle(
             -sum_y,
             sum_x,
             pi
           )
         );
       }

       const double difference_norm = std::hypot(
         difference_x,
         difference_y
       );

       if (difference_norm > zero_threshold) {
         critical_angles.push_back(
           projective_angle(
             -difference_y,
             difference_x,
             pi
           )
         );
       }
     }

     if (i % 256 == 0) {
       Rcpp::checkUserInterrupt();
     }
   }

   const int number_of_raw_angles = static_cast<int>(
     critical_angles.size()
   );

   std::sort(
     critical_angles.begin(),
     critical_angles.end()
   );

   // Remove duplicate interior angles. Angles numerically equal to 0 or pi
   // are represented by the explicit projective boundary at zero.
   std::vector<double> interior_angles;
   interior_angles.reserve(critical_angles.size());

   for (const double angle : critical_angles) {
     if (angle <= tol || angle >= pi - tol) {
       continue;
     }

     if (
         interior_angles.empty() ||
           angle - interior_angles.back() > tol
     ) {
       interior_angles.push_back(angle);
     }
   }

   const int number_of_cells =
     static_cast<int>(interior_angles.size()) + 1;

   Rcpp::NumericVector boundaries(number_of_cells + 1);
   boundaries[0] = 0.0;

   for (int index = 0;
        index < static_cast<int>(interior_angles.size());
        ++index) {
     boundaries[index + 1] = interior_angles[index];
   }

   boundaries[number_of_cells] = pi;

   // Every row contains the 1-based ordering for one angular cell.
   Rcpp::IntegerMatrix orderings(number_of_cells, n);
   Rcpp::LogicalVector has_ties(number_of_cells);

   std::vector<double> squared_projections(n);
   std::vector<int> ordering(n);

   for (int cell = 0; cell < number_of_cells; ++cell) {
     const double alpha_left = boundaries[cell];
     const double alpha_right = boundaries[cell + 1];
     const double theta_mid =
       0.5 * (alpha_left + alpha_right);

     const double direction_x = std::cos(theta_mid);
     const double direction_y = std::sin(theta_mid);

     for (int observation = 0;
          observation < n;
          ++observation) {
       const double projection =
       Y(observation, 0) * direction_x +
       Y(observation, 1) * direction_y;

       const double squared_projection =
         projection * projection;

       if (!std::isfinite(squared_projection)) {
         Rcpp::stop(
           "Squared projections are not finite; the data scale is too large."
         );
       }

       squared_projections[observation] = squared_projection;
       ordering[observation] = observation;
     }

     // Resolve exact numerical ties deterministically by observation index.
     std::sort(
       ordering.begin(),
       ordering.end(),
       [&squared_projections](const int first, const int second) {
         if (
             squared_projections[first] <
               squared_projections[second]
         ) {
           return true;
         }

         if (
             squared_projections[first] >
           squared_projections[second]
         ) {
           return false;
         }

         return first < second;
       }
     );

     double tie_scale = 1.0;

     for (int rank = 0; rank < n; ++rank) {
       const int observation = ordering[rank];

       // R indices start at one.
       orderings(cell, rank) = observation + 1;

       tie_scale = std::max(
         tie_scale,
         std::abs(squared_projections[observation])
       );
     }

     bool cell_has_ties = false;

     for (int rank = 1; rank < n; ++rank) {
       const double previous_value =
         squared_projections[ordering[rank - 1]];

       const double current_value =
         squared_projections[ordering[rank]];

       if (
           current_value - previous_value <=
             tol * tie_scale
       ) {
         cell_has_ties = true;
         break;
       }
     }

     has_ties[cell] = cell_has_ties;

     if (cell % 1024 == 0) {
       Rcpp::checkUserInterrupt();
     }
   }

   return Rcpp::List::create(
     Rcpp::Named("boundaries") = boundaries,
     Rcpp::Named("orderings") = orderings,
     Rcpp::Named("has_ties") = has_ties,
     Rcpp::Named("number_of_raw_angles") = number_of_raw_angles,
     Rcpp::Named("number_of_unique_interior_angles") =
       static_cast<int>(interior_angles.size()),
       Rcpp::Named("number_of_cells") = number_of_cells,
       Rcpp::Named("n") = n
   );
 }
