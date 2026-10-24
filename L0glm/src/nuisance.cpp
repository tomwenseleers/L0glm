// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include <cmath>
#include <limits>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "nuisance.h"
#include "boundary.h"
using Eigen::MatrixXd;using Eigen::VectorXd;
inline double sx_b_loss(const VectorXd& y,const VectorXd& mu,const VectorXd& prior) {
  double value=0;
  for(int i=0;i<y.size();++i) if(prior[i]>0) {
    if(mu[i]<0 || (y[i]>0 && mu[i]<=0) || !std::isfinite(mu[i])) return INFINITY;
    double term=mu[i];
    if(y[i]>0) {
      const double r=(y[i]-mu[i])/mu[i];
      term=std::abs(r)<1e-4 ? .5*mu[i]*r*r*(1-r/3+r*r/6-r*r*r/10+r*r*r*r/15) :
        y[i]*std::log(y[i]/mu[i])-y[i]+mu[i];
    }
    value+=prior[i]*term;
  }
  return value;
}


// Convex nonnegative Poisson background fit with fixed peak means.
// The observed q-dimensional Hessian avoids Fisher stalling at zero counts.
// Only q x q matrices are factorised; channels are independent OpenMP tasks.
// [[Rcpp::export]]
Rcpp::List nuisance_poisson_cpp(const Eigen::Map<Eigen::MatrixXd>& B,const Eigen::Map<Eigen::MatrixXd>& Y,
 const Eigen::Map<Eigen::MatrixXd>& Signal,const Eigen::Map<Eigen::MatrixXd>& Prior,
 const Eigen::Map<Eigen::MatrixXd>& Start,int maxit=60,double tol=1e-8,int nthreads=4) {
  const int n=Y.rows(),m=Y.cols(),q=B.cols();
  if(B.rows()!=n || Signal.rows()!=n || Signal.cols()!=m || Prior.rows()!=n || Prior.cols()!=m ||
     Start.rows()!=q || Start.cols()!=m || maxit<1 || nthreads<1 || tol<=0)
    Rcpp::stop("Incompatible conditional background inputs.");
  if(!B.allFinite() || !Y.allFinite() || !Signal.allFinite() || !Prior.allFinite() || !Start.allFinite() ||
     B.minCoeff()<0 || Y.minCoeff()<0 || Signal.minCoeff()<0 || Prior.minCoeff()<0 || Start.minCoeff()<0)
    Rcpp::stop("Nonnegative finite mean components and weights required.");
  MatrixXd coefficients=Start,background(n,m);VectorXd iterations=VectorXd::Zero(m),kkt=VectorXd::Zero(m);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nthreads) if(nthreads>1 && m>1)
#endif
  for(int channel=0;channel<m;++channel) {
    VectorXd y=Y.col(channel),signal=Signal.col(channel),prior=Prior.col(channel),a=Start.col(channel);
    VectorXd mu=signal+B*a;
    double old=sx_b_loss(y,mu,prior);
    if(!std::isfinite(old)) {
      const double mean=prior.sum()>0 ? prior.dot(y)/prior.sum() : 0;
      a.array()+=std::max(1e-8,mean*.01);mu=signal+B*a;old=sx_b_loss(y,mu,prior);
    }
    for(int iter=0;iter<maxit;++iter) {
      iterations[channel]=iter+1;
      VectorXd w(n),score(n);
      for(int i=0;i<n;++i) {
        const double mean=std::max(mu[i],std::numeric_limits<double>::epsilon());
        w[i]=y[i]>0 ? prior[i]*(y[i]/mean)/mean : 0;score[i]=prior[i]*(1-y[i]/mean);
      }
      MatrixXd BW=B.array().colwise()*w.array();
      MatrixXd G=B.transpose()*BW;
      VectorXd gradient=B.transpose()*score,h=G*a-gradient;
      VectorXd proposal=L0glm::nuisance_qp(G,h,std::vector<char>(q,1)),next=a,pred=mu;
      bool accepted=false;
      for(int attempt=0;attempt<2 && !accepted;++attempt) {
        if(attempt) for(int j=0;j<q;++j)
          proposal[j]=G(j,j)>0 ? std::max(0.0,a[j]-gradient[j]/G(j,j)) : a[j];
        double fraction=1;
        for(int half=0;half<35;++half) {
          next=a+fraction*(proposal-a);pred=signal+B*next;
          const double value=sx_b_loss(y,pred,prior);
          if(std::isfinite(value) && value<=old+1e-12*std::max(1.0,std::abs(old))) {
            accepted=true;old=value;break;
          }
          fraction*=.5;
        }
      }
      if(!accepted) break;
      a=next;mu=pred;
      const VectorXd lower=VectorXd::Zero(q);
      L0glm::boundary_identity(B,y,prior,a,mu,lower,1,true);
      mu=signal+B*a;old=sx_b_loss(y,mu,prior);
      VectorXd derivative(n);
      for(int i=0;i<n;++i) derivative[i]=prior[i]*(1-y[i]/std::max(mu[i],std::numeric_limits<double>::epsilon()));
      const VectorXd gradient_final=B.transpose()*derivative;
      double residual=0;
      for(int j=0;j<q;++j) residual=std::max(residual,
        (a[j]>0 ? std::abs(gradient_final[j]) : std::max(0.0,-gradient_final[j]))/
          std::max(1.0,prior.dot(B.col(j))));
      if(residual<tol) break;
    }
    VectorXd score(n);
    for(int i=0;i<n;++i) {
      const double mean=std::max(mu[i],std::numeric_limits<double>::epsilon());
      score[i]=prior[i]*(1-y[i]/mean);
    }
    VectorXd g=B.transpose()*score;
    double residual=0,scale=1;
    for(int j=0;j<q;++j) {
      const double mass=prior.dot(B.col(j));scale=std::max(scale,mass);
      residual=std::max(residual,a[j]>0 ? std::abs(g[j]) : std::max(0.0,-g[j]));
    }
    coefficients.col(channel)=a;background.col(channel)=B*a;kkt[channel]=residual/scale;
  }
  return Rcpp::List::create(Rcpp::Named("coefficients")=coefficients,Rcpp::Named("background")=background,
    Rcpp::Named("iterations")=iterations,Rcpp::Named("kkt")=kkt);
}

// [[Rcpp::export]]
Eigen::VectorXd nuisance_poisson_deviance_cpp(const Eigen::Map<Eigen::MatrixXd>& Y,const Eigen::Map<Eigen::MatrixXd>& Mu,
 const Eigen::Map<Eigen::MatrixXd>& W,int nthreads=4) {
  if(Y.rows()!=Mu.rows() || Y.cols()!=Mu.cols() || Y.rows()!=W.rows() || Y.cols()!=W.cols())
    Rcpp::stop("Incompatible Poisson scoring dimensions.");
  VectorXd out(Y.cols());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
  for(int j=0;j<Y.cols();++j) out[j]=2*sx_b_loss(Y.col(j),Mu.col(j),W.col(j));
  return out;
}

// [[Rcpp::export]]
Eigen::VectorXd nuisance_quadratic_cpp(const Eigen::Map<Eigen::MatrixXd>& H,
 const Eigen::Map<Eigen::VectorXd>& h,Rcpp::LogicalVector constrained) {
 if(H.rows()!=h.size()||H.cols()!=h.size()||constrained.size()!=h.size()||
    !H.allFinite()||!h.allFinite()) Rcpp::stop("Invalid nuisance quadratic.");
 std::vector<char> bound(h.size());
 for(int j=0;j<h.size();++j) {if(constrained[j]==NA_LOGICAL) Rcpp::stop("Invalid nuisance bounds.");bound[j]=constrained[j];}
 return L0glm::nuisance_qp(H,h,bound);
}

// [[Rcpp::export]]
Rcpp::List nuisance_gaussian_cpp(const Eigen::Map<Eigen::MatrixXd>& B,const Eigen::Map<Eigen::MatrixXd>& Y,
 const Eigen::Map<Eigen::MatrixXd>& Signal,const Eigen::Map<Eigen::MatrixXd>& W,
 const Eigen::Map<Eigen::MatrixXd>& Start,Rcpp::LogicalVector constrained,int nthreads=1) {
 const int n=Y.rows(),m=Y.cols(),q=B.cols();
 if(B.rows()!=n||Signal.rows()!=n||Signal.cols()!=m||W.rows()!=n||W.cols()!=m||
   Start.rows()!=q||Start.cols()!=m||constrained.size()!=q||nthreads<1||
   !B.allFinite()||!Y.allFinite()||!Signal.allFinite()||!W.allFinite()||W.minCoeff()<0)
   Rcpp::stop("Invalid Gaussian nuisance inputs.");
 std::vector<char> bound(q);for(int j=0;j<q;++j) bound[j]=constrained[j];
 MatrixXd coefficients=Start,contribution(n,m);VectorXd kkt=VectorXd::Zero(m);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nthreads) if(nthreads>1 && m>1)
#endif
 for(int c=0;c<m;++c) {
   MatrixXd BW=B.array().colwise()*W.col(c).array();
   MatrixXd H=B.transpose()*BW;VectorXd h=BW.transpose()*(Y.col(c)-Signal.col(c));
   VectorXd a=L0glm::nuisance_qp(H,h,bound),g=H*a-h;
   double residual=0,scale=std::max(1.0,h.lpNorm<Eigen::Infinity>());
   for(int j=0;j<q;++j) residual=std::max(residual,bound[j]&&a[j]==0 ? std::max(0.0,-g[j]) : std::abs(g[j]));
   coefficients.col(c)=a;contribution.col(c)=B*a;kkt[c]=residual/scale;
 }
 return Rcpp::List::create(Rcpp::Named("coefficients")=coefficients,Rcpp::Named("background")=contribution,
    Rcpp::Named("iterations")=Rcpp::IntegerVector(m,1),Rcpp::Named("kkt")=kkt);
}
