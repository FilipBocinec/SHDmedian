test_that("valid bivariate data are accepted", {
  X <- matrix(1:8, ncol = 2L)

  result <- as_bivariate_matrix_d2(X)

  expect_equal(result, X)
  expect_type(result, "double")
  expect_equal(dim(result), c(4L, 2L))
})


test_that("invalid data are rejected", {
  expect_error(
    as_bivariate_matrix_d2(matrix(1:9, ncol = 3L)),
    "exactly two columns"
  )

  expect_error(
    as_bivariate_matrix_d2(matrix(c(1, 2), nrow = 1L)),
    "at least two observations"
  )

  expect_error(
    as_bivariate_matrix_d2(matrix(c(1, 2, 3, Inf), ncol = 2L)),
    "only finite values"
  )

  expect_error(
    as_bivariate_matrix_d2(matrix(c("a", "b", "c", "d"), ncol = 2L)),
    "numeric matrix"
  )
})


test_that("center validation works", {
  expect_equal(
    as_center_vector_d2(c(1, 2)),
    c(1, 2)
  )

  expect_error(
    as_center_vector_d2(1),
    "length two"
  )

  expect_error(
    as_center_vector_d2(c(1, NA_real_)),
    "finite numeric vector"
  )
})
