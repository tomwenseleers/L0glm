#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include <L0glm/scalar_families.h>
#include <L0glm/nuisance.h>
#include <set>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

// Fisher scoring on a linear-interpolated concave log profile. A small dual
// nonnegative QP enforces the slope inequalities; likelihood is original scale.
// [[Rcpp::export]]
List slcn_log_profile(NumericMatrix y,NumericMatrix mean,NumericVector profile,
 NumericVector spectrum,int family,NumericVector phi,NumericVector theta,
 NumericVector weights,int full_n,int first=0,int maxit=4,double censor=.5,
 int q=128,int nthreads=1,int apex_lo=-1,int apex_hi=-1,int mode=2) {
  const int n=y.nrow(),m=y.ncol();q=std::min(n,q);
  const bool smooth=mode==0,whole=mode!=1;
  if(n<3||q<3||first<0||maxit<0||nthreads<1||family<0||family>6||mode<0||mode>2||
     mean.nrow()!=n||mean.ncol()!=m||profile.size()!=n||spectrum.size()!=m||
     full_n<first+n||!(weights.size()==1||weights.size()==full_n*m)||
     !(phi.size()==1||phi.size()==m)||!(theta.size()==1||theta.size()==m))
    stop("Invalid log-profile inputs.");
  for(double a:profile)if(!std::isfinite(a)||a<=0)stop("Positive finite profile required.");
  auto w=[&](int i,int j){return weights[weights.size()==1?0:first+i+full_n*j];};
  auto p=[&](int j){return phi[phi.size()==1?0:j];};
  auto th=[&](int j){return theta[theta.size()==1?0:j];};
  std::vector<int> knots;Eigen::MatrixXd C,D;bool lo=false,hi=false;
  if(smooth) {
    // Integrating nonnegative, piecewise-constant curvature twice gives a
    // smooth concave log profile. The intercept and slope remain unrestricted.
    C.resize(n,q);const int J=q-2;
    for(int i=0;i<n;++i) {
      const double t=(double)i/(n-1);C(i,0)=1;C(i,1)=t;
      for(int j=0;j<J;++j) {
        const double a=std::max(0.,t-(double)j/J),b=std::max(0.,t-(double)(j+1)/J);
        C(i,j+2)=-.5*(a*a-b*b);
      }
    }
    for(int j=0;j<q;++j)C.col(j)/=C.col(j).norm();
  } else {
  std::set<int> index;const int regular=std::max(3,q/3);
  for(int j=0;j<regular;++j) index.insert((int)std::round((double)j*(n-1)/(regular-1)));
  {
    double total=0;for(double a:profile)total+=a;
    int i=0;double cumulative=profile[0];
    const int by_mass=q/3;
    for(int j=1;j<by_mass;++j) {
      const double target=total*j/by_mass;
      while(i<n-1&&cumulative<target)cumulative+=profile[++i];
      index.insert(i);
    }
    // Include current sharp curvature points so an existing kink can relax.
    std::vector<std::pair<double,int>> curvature;double peak=0;
    for(int j=1;j<n-1;++j) {
      double a=std::max(0.,-std::log(profile[j-1])+2*std::log(profile[j])-std::log(profile[j+1]));
      curvature.push_back({a,j});peak=std::max(peak,a);
    }
    std::sort(curvature.begin(),curvature.end(),[](const auto& a,const auto& b){return a.first>b.first;});
    int significant=0;for(const auto& a:curvature)if(a.first>std::max(1e-9,peak*1e-3))++significant;
    const int limit=significant<=q/2?significant:q/3;
    for(int j=0;j<limit;++j)index.insert(curvature[j].second);
  }
  lo=apex_lo>0&&apex_lo<n;hi=apex_hi>=0&&apex_hi<n-1;
  if(lo){index.insert(apex_lo-1);index.insert(apex_lo);}
  if(hi){index.insert(apex_hi);index.insert(apex_hi+1);}
  knots.assign(index.begin(),index.end());q=knots.size();
  C=Eigen::MatrixXd::Zero(n,q);D=Eigen::MatrixXd::Zero(q-2+lo+hi,q);
  int left=0;
  for(int i=0;i<n;++i) {
    while(left<q-2&&i>knots[left+1])++left;
    const double fraction=(double)(i-knots[left])/(knots[left+1]-knots[left]);
    C(i,left)=1-fraction;C(i,left+1)=fraction;
  }
  for(int j=0;j<q-2;++j) {
    const double a=1.0/(knots[j+1]-knots[j]),b=1.0/(knots[j+2]-knots[j+1]);
    D(j,j)=a;D(j,j+1)=-a-b;D(j,j+2)=b;D.row(j)/=D.row(j).norm();
  }
  int next_constraint=q-2;
  if(lo){D.row(next_constraint)=C.row(apex_lo-1)-C.row(apex_lo);D.row(next_constraint)/=D.row(next_constraint).norm();++next_constraint;}
  if(hi){D.row(next_constraint)=C.row(apex_hi+1)-C.row(apex_hi);D.row(next_constraint)/=D.row(next_constraint).norm();}
  }
  NumericMatrix rest(n,m);NumericVector h=clone(profile);
  Eigen::VectorXd eta(n),gradient(n),information(n),target(n),W(n);
  for(int i=0;i<n;++i) {
    if(!std::isfinite(h[i])||h[i]<=0)stop("Positive finite profile required.");
    eta[i]=std::log(h[i]);
  }
  for(int j=0;j<m;++j)for(int i=0;i<n;++i)rest(i,j)=std::max(0.,mean(i,j)-profile[i]*spectrum[j]);
  auto change=[&](const NumericVector& old,const NumericVector& next) {
    std::vector<double> sums(n);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
    for(int i=0;i<n;++i) {
      long double sum=0;
      for(int j=0;j<m;++j)if(w(i,j)>0&&spectrum[j]>0)
        sum+=w(i,j)*L0glmScalar::change(y(i,j),rest(i,j)+old[i]*spectrum[j],
          rest(i,j)+next[i]*spectrum[j],family,p(j),th(j),censor);
      sums[i]=(double)sum;
    }
    long double sum=0;for(double a:sums)sum+=a;return (double)sum;
  };
  std::vector<double> history;double residual=0;
  for(int iteration=0;iteration<maxit;++iteration) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
    for(int i=0;i<n;++i) {
      long double g=0,I=0;
      for(int j=0;j<m;++j)if(w(i,j)>0&&spectrum[j]>0) {
        const double signal=h[i]*spectrum[j],mu=std::max(1e-15,rest(i,j)+signal);
        const auto a=L0glmScalar::derivatives(y(i,j),mu,family,p(j),th(j),censor);
        g+=w(i,j)*signal*a.gradient;I+=w(i,j)*signal*signal*a.information;
      }
      gradient[i]=(double)g;information[i]=(double)I;
    }
    const double top=information.maxCoeff();if(!(top>0))break;
    for(int i=0;i<n;++i) {
      W[i]=information[i]/top+1e-8;
      target[i]=(whole?eta[i]:0)-gradient[i]/(smooth?std::max(information[i],top*1e-8):(information[i]+top*1e-8));
    }
    Eigen::MatrixXd G=C.transpose()*W.asDiagonal()*C;
    G.diagonal().array()+=std::max(1e-12,1e-10*G.diagonal().maxCoeff());
    Eigen::VectorXd proposed;
    if(smooth) {
      Eigen::VectorXd rhs=C.transpose()*W.asDiagonal()*target;
      std::vector<char> bounds(q,1);bounds[0]=bounds[1]=0;
      proposed=C*L0glm::nuisance_qp(G,rhs,bounds);
    } else {
      Eigen::VectorXd rhs=C.transpose()*(W.array()*target.array()).matrix();
      Eigen::LLT<Eigen::MatrixXd> factor(G);
      if(factor.info()!=Eigen::Success)stop("Log-profile Gram failed.");
      Eigen::VectorXd unconstrained=factor.solve(rhs);
      Eigen::MatrixXd inverse_D=factor.solve(D.transpose());
      Eigen::MatrixXd dual_G=D*inverse_D;dual_G=(dual_G+dual_G.transpose()).eval()*.5;
      Eigen::VectorXd bound=Eigen::VectorXd::Zero(D.rows());
      if(!whole)for(int j=0;j<q-2;++j) {
        const int k=knots[j+1];
        const double a=1.0/(knots[j+1]-knots[j]),b=1.0/(knots[j+2]-knots[j+1]);
        bound[j]=std::max(0.,-(eta[k-1]-2*eta[k]+eta[k+1]))/std::sqrt(a*a+(a+b)*(a+b)+b*b);
      }
      if(lo)bound[q-2]=-1e-6;
      if(hi)bound[q-2+lo]=-1e-6;
      Eigen::VectorXd dual_rhs=D*unconstrained-bound;
      Eigen::VectorXd multiplier=L0glm::nuisance_qp(dual_G,dual_rhs,std::vector<char>(D.rows(),1));
      Eigen::VectorXd node=unconstrained-inverse_D*multiplier;
      residual=std::max(0.,(D*node-bound).maxCoeff());
      proposed=C*node;
      if(!whole)proposed+=eta;
    }
    // Remove numerical inequality residuals on the full grid.
    if(!smooth) {
    std::vector<int> hull;
    for(int i=0;i<n;++i) {
      hull.push_back(i);
      while(hull.size()>2) {
        int k=hull.size(),a=hull[k-3],b=hull[k-2],c=hull[k-1];
        if((proposed[b]-proposed[a])/(b-a)>=
           (proposed[c]-proposed[b])/(c-b))break;
        hull[k-2]=c;hull.pop_back();
      }
    }
    Eigen::VectorXd feasible=proposed;
    for(size_t j=1;j<hull.size();++j)for(int i=hull[j-1];i<=hull[j];++i) {
      const int a=hull[j-1],b=hull[j];
      feasible[i]=proposed[a]+(proposed[b]-proposed[a])*(i-a)/(b-a);
    }
    proposed=feasible;
    }
    bool accepted=false;NumericVector next(n);double delta=0;
    for(int half=0;half<=(smooth?20:22);++half) {
      bool valid=true;const double fraction=std::ldexp(1.,-half);
      for(int i=0;i<n;++i) {
        const double value=eta[i]+fraction*(proposed[i]-eta[i]);
        if(value>600||value< -600||!std::isfinite(value)){valid=false;break;}
        next[i]=std::exp(value);
      }
      if(!valid)continue;
      delta=change(h,next);
      if(std::isfinite(delta)&&delta<=1e-7){accepted=true;break;}
    }
    if(!accepted)break;
    h=next;for(int i=0;i<n;++i)eta[i]=std::log(h[i]);history.push_back(delta);
    if(-delta<1e-6)break;
    checkUserInterrupt();
  }
  return List::create(_["profile"]=h,_["history"]=history,_["knots"]=knots,_["dual_residual"]=residual);
}
