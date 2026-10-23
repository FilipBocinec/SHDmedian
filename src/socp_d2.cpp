// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>
#include <numeric>

namespace point_socp {
const double PI = std::acos(-1.0);
struct Vec { double x, y; };
struct Event { double angle; int j; };
struct Group { double angle; std::vector<int> ids; };
struct PointGeometry {
  std::vector<double> angle;
  std::vector<int> less;
  int equal;
};
struct Geometry {
  std::vector<Vec> Y;
  std::vector<PointGeometry> points;
  bool cached;
  double angle_tol;
};
struct Arc { int i; double left, right; };
struct Cut { double angle, rhs; };
struct Block { int i; std::vector<Vec> normal; std::vector<double> offset; };
struct Level {
  std::vector<Vec> Y;
  std::vector<Arc> lower, upper;
  std::vector<Block> blocks;
  std::vector<Cut> cuts;
  bool structural_empty = false;
  int k;
};

double canonical(double t, double tol=0) {
  t = std::fmod(t, PI);
  if (t < 0) t += PI;
  if (t < tol || PI-t < tol) t = 0;
  return t;
}
double projection2(const Vec &v, double c, double s) {
  double a = v.x*c+v.y*s; return a*a;
}
bool same_outer(const Vec &a, const Vec &b) {
  return (a.x == b.x && a.y == b.y) || (a.x == -b.x && a.y == -b.y);
}
std::vector<Group> groups(std::vector<Event> &events, double tol) {
  for (Event &e: events) e.angle = canonical(e.angle,tol);
  std::sort(events.begin(),events.end(),[](const Event &a,const Event &b){
    return a.angle < b.angle;
  });
  std::vector<Group> out;
  for (const Event &e: events) {
    if (out.empty() || e.angle-out.back().angle > tol)
      out.push_back(Group{e.angle,{}});
    if (e.j >= 0) out.back().ids.push_back(e.j);
  }
  if (out.empty()) out.push_back(Group{0,{}});
  for (Group &g: out) {
    std::sort(g.ids.begin(),g.ids.end());
    g.ids.erase(std::unique(g.ids.begin(),g.ids.end()),g.ids.end());
  }
  return out;
}
double midpoint(const std::vector<Group> &g, int t) {
  return (g[t].angle + (t+1<(int)g.size()?g[t+1].angle:g[0].angle+PI))/2;
}
PointGeometry point_geometry(const Geometry &G, int i) {
  int n=G.Y.size(), equal=0;
  const Vec &yi=G.Y[i];
  std::vector<Event> events; events.reserve(2*n);
  for (int j=0;j<n;++j) {
    const Vec &yj=G.Y[j];
    if (same_outer(yi,yj)) {++equal; continue;}
    Vec minus{yi.x-yj.x,yi.y-yj.y}, plus{yi.x+yj.x,yi.y+yj.y};
    events.push_back(Event{std::atan2(minus.x,-minus.y),j});
    events.push_back(Event{std::atan2(plus.x,-plus.y),j});
  }
  std::vector<Group> g=groups(events,G.angle_tol);
  PointGeometry ans; ans.equal=equal;
  ans.angle.reserve(g.size()); ans.less.reserve(g.size());
  std::vector<unsigned char> below(n,0);
  double t=midpoint(g,0), c=std::cos(t), s=std::sin(t);
  double ri=projection2(yi,c,s); int less=0;
  for(int j=0;j<n;++j) {
    if(!same_outer(yi,G.Y[j])) below[j]=projection2(G.Y[j],c,s)<ri;
    less+=below[j];
  }
  for(int h=0;h<(int)g.size();++h) {
    if(h>0) {
      t=midpoint(g,h); c=std::cos(t); s=std::sin(t); ri=projection2(yi,c,s);
      for(int j:g[h].ids) {
        bool b=projection2(G.Y[j],c,s)<ri;
        less+=(int)b-(int)below[j]; below[j]=b;
      }
    }
    ans.angle.push_back(g[h].angle); ans.less.push_back(less);
  }
  return ans;
}
std::vector<Arc> runs(const PointGeometry &p, const std::vector<unsigned char> &flag,
                      int i, bool value) {
  int m=flag.size(); std::vector<Arc> ans;
  int count=0; for(bool f:flag) count+=(f==value);
  if(count==0) return ans;
  if(count==m) {ans.push_back(Arc{i,0,PI}); return ans;}
  for(int j=0;j<m;++j) if(flag[j]==value && flag[(j+m-1)%m]!=value) {
    int end=(j+1)%m;
    while(flag[end]==value) end=(end+1)%m;
    double a=p.angle[j], b=p.angle[end]; if(end<=j) b+=PI;
    ans.push_back(Arc{i,a,b});
  }
  return ans;
}
void add_cut(Level &L,int i,double theta,double scale) {
  theta=canonical(theta);
  L.cuts.push_back(Cut{theta,projection2(L.Y[i],std::cos(theta),std::sin(theta))/scale});
}
Level* build_level(const Geometry &G, int k, double scale) {
  Level *L=new Level; L->Y=G.Y; L->k=k;
  int n=G.Y.size();
  for(int i=0;i<n;++i) {
    if((i&31)==0) Rcpp::checkUserInterrupt();
    PointGeometry tmp;
    const PointGeometry *p;
    if(G.cached) p=&G.points[i]; else {tmp=point_geometry(G,i); p=&tmp;}
    int m=p->less.size();
    std::vector<unsigned char> low(m),up(m);
    for(int j=0;j<m;++j) {
      low[j]=p->less[j]<=k-1;
      up[j]=n-p->equal-p->less[j]<=k-1;
    }
    auto lower=runs(*p,low,i,true), gaps=runs(*p,low,i,false);
    auto upper=runs(*p,up,i,true);
    L->lower.insert(L->lower.end(),lower.begin(),lower.end());
    L->upper.insert(L->upper.end(),upper.begin(),upper.end());
    bool zero=(G.Y[i].x==0 && G.Y[i].y==0);
    if(!lower.empty() && !zero) {
      Block block; block.i=i;
      for(const Arc &gap:gaps) {
        block.normal.push_back(Vec{std::cos(gap.left+gap.right),std::sin(gap.left+gap.right)});
        block.offset.push_back(std::cos(gap.right-gap.left));
      }
      L->blocks.push_back(std::move(block));
    }
    for(const Arc &arc:upper) {
      // On each side of Yi^perp, |u'Yi| is linear on the cone generated by
      // endpoint rays. Sublinearity of sqrt(u'Sigma u), Sigma >= 0, makes
      // the two endpoint inequalities sufficient for the whole arc.
      if(zero || arc.right-arc.left >= PI-G.angle_tol) {
        L->structural_empty=true; continue;
      }
      double perp=canonical(std::atan2(G.Y[i].x,-G.Y[i].y));
      double distance=canonical(perp-canonical(arc.left));
      if(distance<=arc.right-arc.left+G.angle_tol) {
        L->structural_empty=true; continue;
      }
      add_cut(*L,i,arc.left,scale); add_cut(*L,i,arc.right,scale);
    }
  }
  // Keep nearby-but-distinct cuts. Drop only exact duplicate angles,
  // retaining the tighter right hand side.
  std::sort(L->cuts.begin(),L->cuts.end(),[](const Cut&a,const Cut&b){return a.angle<b.angle;});
  std::vector<Cut> distinct;
  for(const Cut &cut:L->cuts) {
    if(!distinct.empty() && cut.angle==distinct.back().angle)
      distinct.back().rhs=std::min(distinct.back().rhs,cut.rhs);
    else distinct.push_back(cut);
  }
  L->cuts.swap(distinct); return L;
}
double quadratic(double a,double c,double b,double theta) {
  double u=std::cos(theta),v=std::sin(theta);
  return a*u*u+2*c*u*v+b*v*v;
}
std::pair<double,double> arc_min(double a,double c,double b,double left,double right) {
  double val=quadratic(a,c,b,left), arg=left;
  double v=quadratic(a,c,b,right); if(v<val){val=v;arg=right;}
  double t=canonical((std::atan2(2*c,a-b)+PI)/2);
  double shifted=left+canonical(t-canonical(left));
  if(shifted<=right+1e-14) {
    v=(a+b)/2-std::hypot((a-b)/2,c);
    if(v<val){val=v;arg=shifted;}
  }
  return std::make_pair(val,canonical(arg));
}
}

// [[Rcpp::export]]
Rcpp::List shd_socp_prepare_cpp(Rcpp::NumericMatrix Y, bool cache=true,
                              double angle_tol=1e-12) {
  using namespace point_socp;
  if(Y.ncol()!=2 || Y.nrow()<1) Rcpp::stop("Y must have two columns.");
  Rcpp::XPtr<Geometry> ptr(new Geometry,true);
  ptr->cached=cache; ptr->angle_tol=angle_tol;
  for(int i=0;i<Y.nrow();++i) ptr->Y.push_back(Vec{Y(i,0),Y(i,1)});
  double bytes=0, events=0;
  if(cache) {
    ptr->points.reserve(Y.nrow());
    for(int i=0;i<Y.nrow();++i) {
      if((i&31)==0) Rcpp::checkUserInterrupt();
      ptr->points.push_back(point_geometry(*ptr,i));
      events+=ptr->points.back().angle.size();
      bytes+=ptr->points.back().angle.size()*(sizeof(double)+sizeof(int));
    }
  }
  return Rcpp::List::create(Rcpp::_["ptr"]=ptr,Rcpp::_["cache_bytes"]=bytes,
                            Rcpp::_["event_count"]=events);
}

// [[Rcpp::export]]
Rcpp::List shd_socp_model_cpp(SEXP geometry_ptr, int k, double scale=1) {
  using namespace point_socp;
  Rcpp::XPtr<Geometry> G(geometry_ptr);
  if(k<1 || k>(int)G->Y.size() || !(scale>0)) Rcpp::stop("Invalid k or scale.");
  Rcpp::XPtr<Level> L(build_level(*G,k,scale),true);
  int ngap=0; for(const Block &b:L->blocks) ngap+=b.offset.size();
  int nupper=L->cuts.size(), nblock=L->blocks.size();
  Rcpp::List info=Rcpp::List::create(
    Rcpp::_["level_ptr"]=L,Rcpp::_["structural_empty"]=L->structural_empty,
    Rcpp::_["n_lower_soc"]=nblock,Rcpp::_["n_multipliers"]=ngap,
    Rcpp::_["n_upper_cuts"]=nupper,Rcpp::_["n_lower_arcs"]=(int)L->lower.size(),
    Rcpp::_["n_upper_arcs"]=(int)L->upper.size());
  if(L->structural_empty) return info;
  int nvar=4+ngap, nlin=1+ngap+nupper, nrow=nlin+3*(1+nblock);
  std::vector<int> ii,jj; std::vector<double> xx;
  ii.reserve(10*nblock+5*ngap+3*nupper+10); jj.reserve(ii.capacity());xx.reserve(ii.capacity());
  Rcpp::NumericVector rhs(nrow),obj(nvar); obj[3]=-1;
  auto add=[&](int row,int col,double value){
    if(value!=0){ii.push_back(row+1);jj.push_back(col+1);xx.push_back(value);}
  };
  add(0,3,-1); // delta >= 0
  for(int j=0;j<ngap;++j) add(1+j,4+j,-1);
  for(int j=0;j<nupper;++j) {
    double c=std::cos(L->cuts[j].angle),s=std::sin(L->cuts[j].angle);
    int row=1+ngap+j;
    add(row,0,c*c); add(row,1,2*c*s); add(row,2,s*s);
    rhs[row]=L->cuts[j].rhs;
  }
  int row=nlin;
  add(row,0,-1);add(row,2,-1);add(row,3,2);
  add(row+1,0,-1);add(row+1,2,1);add(row+2,1,-2);
  row+=3; int first=4;
  for(const Block &block:L->blocks) {
    Vec yi=L->Y[block.i];
    double a=yi.x*yi.x/scale,c=yi.x*yi.y/scale,b=yi.y*yi.y/scale;
    rhs[row]=-(a+b)/2;rhs[row+1]=-(a-b)/2;rhs[row+2]=-c;
    add(row,0,-.5);add(row,2,-.5);
    add(row+1,0,-.5);add(row+1,2,.5);add(row+2,1,-1);
    for(int j=0;j<(int)block.offset.size();++j) {
      add(row,first+j,block.offset[j]);
      add(row+1,first+j,-block.normal[j].x);
      add(row+2,first+j,-block.normal[j].y);
    }
    first+=block.offset.size();row+=3;
  }
  info["i"]=Rcpp::wrap(ii); info["j"]=Rcpp::wrap(jj);info["x"]=Rcpp::wrap(xx);
  info["b"]=rhs;info["q"]=obj;info["nrow"]=nrow;info["nvar"]=nvar;
  info["nlin"]=nlin;info["nsoc"]=nblock+1;
  return info;
}

// [[Rcpp::export]]
Rcpp::List shd_socp_validate_cpp(SEXP level_ptr, Rcpp::NumericMatrix Sigma) {
  using namespace point_socp;
  Rcpp::XPtr<Level> L(level_ptr);
  double lo=std::numeric_limits<double>::infinity(),hi=lo,lo_t=0,hi_t=0;
  for(int side=0;side<2;++side) {
    const auto &arcs=(side==0?L->lower:L->upper);
    for(const Arc &arc:arcs) {
      Vec v=L->Y[arc.i]; double sign=(side==0?1:-1);
      auto ans=arc_min(sign*(Sigma(0,0)-v.x*v.x),sign*(Sigma(0,1)-v.x*v.y),
                       sign*(Sigma(1,1)-v.y*v.y),arc.left,arc.right);
      if(side==0 && ans.first<lo){lo=ans.first;lo_t=ans.second;}
      if(side==1 && ans.first<hi){hi=ans.first;hi_t=ans.second;}
    }
  }
  return Rcpp::List::create(Rcpp::_["lower_slack"]=lo,Rcpp::_["upper_slack"]=hi,
    Rcpp::_["lower_angle"]=lo_t,Rcpp::_["upper_angle"]=hi_t);
}

// [[Rcpp::export]]
Rcpp::List shd_depth_point_cpp(Rcpp::NumericMatrix Y,Rcpp::NumericMatrix Sigma,
                             double tolerance=0,double angle_tol=1e-13) {
  using namespace point_socp;
  int n=Y.nrow();
  std::vector<Vec> beta(n);std::vector<double> alpha(n);
  std::vector<Event> events;events.reserve(4*n);
  for(int i=0;i<n;++i) {
    double a=Sigma(0,0)-Y(i,0)*Y(i,0),b=Sigma(1,1)-Y(i,1)*Y(i,1);
    alpha[i]=(a+b)/2;beta[i]=Vec{(a-b)/2,Sigma(0,1)-Y(i,0)*Y(i,1)};
    double radius=std::hypot(beta[i].x,beta[i].y);
    if(radius==0) continue;
    double phi=std::atan2(beta[i].y,beta[i].x);
    for(int sign: {-1,1}) {
      double q=(sign*tolerance-alpha[i])/radius;
      if(q>-1 && q<1) {
        double a0=std::acos(q);
        events.push_back(Event{(phi+a0)/2,i});
        events.push_back(Event{(phi-a0)/2,i});
      }
    }
  }
  auto g=groups(events,angle_tol);
  std::vector<unsigned char> inside(n),outside(n);
  auto evaluate=[&](int i,double theta){
    return alpha[i]+beta[i].x*std::cos(2*theta)+beta[i].y*std::sin(2*theta);
  };
  int ni=0,no=0;
  for(int i=0;i<n;++i) {
    double f=evaluate(i,midpoint(g,0));
    inside[i]=f>=-tolerance;outside[i]=f<=tolerance;
    ni+=inside[i];no+=outside[i];
  }
  int best=std::min(ni,no), minin=ni,minout=no;double worst=midpoint(g,0);
  for(int h=1;h<(int)g.size();++h) {
    double t=midpoint(g,h);
    for(int i:g[h].ids) {
      double f=evaluate(i,t);bool in=f>=-tolerance,out=f<=tolerance;
      ni+=(int)in-(int)inside[i];no+=(int)out-(int)outside[i];
      inside[i]=in;outside[i]=out;
    }
    minin=std::min(minin,ni);minout=std::min(minout,no);
    if(std::min(ni,no)<best){best=std::min(ni,no);worst=t;}
  }
  return Rcpp::List::create(Rcpp::_["k"]=best,Rcpp::_["depth"]=(double)best/n,
    Rcpp::_["minimum_inside"]=minin,Rcpp::_["minimum_outside"]=minout,
    Rcpp::_["worst_angle"]=canonical(worst),Rcpp::_["n_events"]=(int)g.size(),
    Rcpp::_["tolerance"]=tolerance);
}
