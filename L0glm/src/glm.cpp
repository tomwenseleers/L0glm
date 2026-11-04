// Channel-level OpenMP owns parallelism; Eigen kernels run within each channel.
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <memory>
#include <cfloat>
#include <type_traits>
#include <L0glm/nuisance.h>
#include "boundary.h"
#include <L0glm/scalar_families.h>
#include "ridge_profile.h"
#include "solver_native.h"
#include <algorithm>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using Eigen::MatrixXd;using Eigen::VectorXd;using Eigen::SparseMatrix;
namespace L0glmNative {
const double eps=std::numeric_limits<double>::epsilon();
struct Family {
 int kind,link,scalar;double theta,parameter,censor,probit_limit,cauchy_limit;
 bool eta_valid(double e) const {
  return std::isfinite(e)&&(scalar>=0||((link!=2&&link!=8)||e>0))&&(scalar>=0||link!=3||e!=0);
 }
 double mean(double e) const {
  switch(link) {
   case 0:return e;case 1:return std::max(eps,std::exp(e));case 2:return e*e;case 3:return 1/e;
   case 4:{double t=e< -30?eps:e>30?1/eps:std::exp(e);return t/(1+t);}
   case 5:return R::pnorm5(std::max(-probit_limit,std::min(probit_limit,e)),0,1,1,0);
   case 6:return std::max(eps,std::min(1-eps,-std::expm1(-std::exp(e))));
   case 7:return R::pcauchy(std::max(-cauchy_limit,std::min(cauchy_limit,e)),0,1,1,0);
   case 8:return 1/std::sqrt(e);
  }
  return NA_REAL;
 }
 double derivative(double e) const {
  switch(link) {
   case 0:return 1;case 1:return std::max(eps,std::exp(e));case 2:return 2*e;case 3:return -1/(e*e);
   case 4:{if(e>30||e< -30)return eps;double t=std::exp(e);return t/((1+t)*(1+t));}
   case 5:return std::max(eps,R::dnorm4(e,0,1,0));
   case 6:{double t=std::exp(std::min(700.,e));return std::max(eps,t*std::exp(-t));}
   case 7:return std::max(eps,R::dcauchy(e,0,1,0));
   case 8:return -1/(2*std::pow(e,1.5));
  }
  return NA_REAL;
 }
 double variance(double mu) const {
  if(scalar>=0)return scalar==4?parameter*mu*mu:parameter*mu;
  switch(kind){case 0:return 1;case 1:return mu;case 2:return mu*(1-mu);
    case 3:return mu*mu;case 4:return mu*mu*mu;case 5:return mu+mu*mu/theta;}
  return NA_REAL;
 }
 double loss(double y,double mu,double eta) const {
  if(scalar>=0) {
   if(mu<0||(mu==0&&!(scalar==5||censor>0)))return R_PosInf;
   return 2*L0glmScalar::terms(y,mu,scalar,parameter,theta,censor).loss;
  }
  if(kind==0)return (y-mu)*(y-mu);
  if(kind==2) {
   if(!(mu>0&&mu<1))return R_PosInf;
   return 2*((y>0?y*std::log(y/mu):0)+(y<1?(1-y)*std::log((1-y)/(1-mu)):0));
  }
  if(mu<0)return R_PosInf;mu=std::max(eps,mu);
  if(kind==1)return 2*((y==0?0:y*std::log(y/mu))-(y-mu));
  if(kind==3)return 2*((y-mu)/mu-std::log(y/mu));
  if(kind==4) {
   if(link==1){double a=(y*std::exp(-eta)-1)/std::sqrt(y);return a*a;}
   return (y-mu)*(y-mu)/(y*mu*mu);
  }
  if(kind==5) {
   if(!(mu>0))return R_PosInf;
   return 2*((y==0?0:y*std::log(y/mu))-(y+theta)*std::log((y+theta)/(mu+theta)));
  }
  return R_PosInf;
 }
 double deviance(const VectorXd& y,const VectorXd& eta,const VectorXd& w) const {
  long double total=0;
  for(int i=0;i<y.size();++i)if(w[i]>0) {
   if(!eta_valid(eta[i]))return R_PosInf;
   double mu=mean(eta[i]);if(!std::isfinite(mu))return R_PosInf;
   double value=loss(y[i],mu,eta[i]);if(!std::isfinite(value))return R_PosInf;total+=w[i]*value;
  }
  return (double)total;
 }
};
struct State {
  VectorXd y,w,offset,beta,eta,mu,previous_eta,anchor,prior,scale,lower,upper,W,z,Wret,zret,lambdas,old_prediction;
 std::vector<char> included;double phi,lambda,tol,thresh,previous_dev;
 int maxit,offset_df,free_intercept,npositive,iterations=0,updates=0,restarts=0,stable=0;
  MatrixXd cacheG;VectorXd cacheXY,cacheW;
  bool estimate,normalized,hard,frozen,cache=false,cache_weights=false,
       converged=false,failed=false,solver_failed=false;
 bool auto_beta=false,auto_phi=false,gaussian_name=false,binomial_name=false,signal_family=false;
 std::vector<char> allowed;double reference_p=1;
 VectorXd center_shift,original_prior;
 int center_intercept=-1;
 double selection_lambda=0;
 bool identity_boundary=false;
 int solver_inexact=0,solver_fallbacks=0;
};
VectorXd predict(const Eigen::Map<MatrixXd>& X,const VectorXd& b,const VectorXd& offset) {
 const int n=X.rows();int present=(b.array()!=0).count();
 if(present>X.cols()/2)return offset+X*b;
 VectorXd value=VectorXd::Zero(n);int inc=1;
 for(int j=0;j<X.cols();++j)if(b[j]!=0)
  F77_CALL(daxpy)(&n,&b[j],X.col(j).data(),&inc,value.data(),&inc);
 return offset+value;
}
VectorXd predict(const Eigen::MappedSparseMatrix<double>& X,const VectorXd& b,const VectorXd& offset) {
 VectorXd value=VectorXd::Zero(X.rows());
 for(int j=0;j<X.cols();++j)if(b[j]!=0)
  for(Eigen::MappedSparseMatrix<double>::InnerIterator it(X,j);it;++it)value[it.row()]+=it.value()*b[j];
 return offset+value;
}
#include "glm_initialize.h"
VectorXd ridge(const Eigen::Map<MatrixXd>& X,const VectorXd& z,const VectorXd& W,
 const std::vector<int>& active,const VectorXd& penalties,const VectorXd& prior,
 const VectorXd& scale,VectorXd previous,bool& failed,MatrixXd& cacheG,VectorXd& cacheXY,
 VectorXd& cacheW,bool cache,bool cache_weights) {
 std::vector<int> good;for(int i=0;i<W.size();++i)if(W[i]>0&&std::isfinite(z[i]))good.push_back(i);
 const int n=good.size(),p=active.size();if(!n){failed=true;return previous;}
  MatrixXd A;VectorXd response(n),root(n),d(p),center(p);
  for(int i=0;i<n;++i){root[i]=std::sqrt(W[good[i]]);response[i]=root[i]*z[good[i]];}
 for(int j=0;j<p;++j){int a=active[j];d[j]=penalties[a]/(scale[a]*scale[a]);center[j]=prior[a]*scale[a];}
  const bool fixed=cache_weights&&X.cols()<=n&&X.cols()<=2000;
  if(fixed&&cacheW.size()&&
      (W-cacheW).cwiseAbs().maxCoeff()>1e-12*std::max(1.,W.cwiseAbs().maxCoeff()))cacheG.resize(0,0);
  const bool use_cache=(cache||fixed)&&p<=n;
  if(!use_cache) {
   A.resize(n,p);for(int j=0;j<p;++j)for(int i=0;i<n;++i)A(i,j)=X(good[i],active[j])*root[i]/scale[active[j]];
 }else if(!cacheG.rows()) {
   int full=X.cols();MatrixXd all(n,full);for(int j=0;j<full;++j)for(int i=0;i<n;++i)all(i,j)=X(good[i],j)*root[i]/scale[j];
  cacheG=MatrixXd::Zero(full,full);cacheXY=VectorXd::Zero(full);double one=1,zero=0;int inc=1;
  F77_CALL(dsyrk)("L","T",&full,&n,&one,all.data(),&n,&zero,cacheG.data(),&full FCONE FCONE);
  F77_CALL(dgemv)("T",&n,&full,&one,all.data(),&n,response.data(),&inc,&zero,cacheXY.data(),&inc FCONE);
   if(fixed)cacheW=W;
 }
 VectorXd fit(p);
 if(p>n) {
   if(!L0glm::ridge_dual(A,response,d,center,fit,false)){failed=true;return previous;}
 }else {
   MatrixXd G=MatrixXd::Zero(p,p);VectorXd rhs(p);double one=1,zero=0;int inc=1;
    if(use_cache) {
    for(int j=0;j<p;++j){rhs[j]=cacheXY[active[j]];for(int k=0;k<=j;++k)
      G(j,k)=cacheG(std::max(active[j],active[k]),std::min(active[j],active[k]));}
     if(fixed&&!cache) {
      const int total=X.rows();int inc=1;VectorXd wz(total);
      for(int i=0;i<total;++i)wz[i]=W[i]>0?W[i]*z[i]:0;
      for(int j=0;j<p;++j)
       rhs[j]=F77_CALL(ddot)(&total,X.col(active[j]).data(),&inc,wz.data(),&inc)/scale[active[j]];
     }
   }else {
    F77_CALL(dsyrk)("L","T",&p,&n,&one,A.data(),&n,&zero,G.data(),&p FCONE FCONE);
    F77_CALL(dgemv)("T",&n,&p,&one,A.data(),&n,response.data(),&inc,&zero,rhs.data(),&inc FCONE);
   }
   for(int j=0;j<p;++j){G(j,j)+=d[j];rhs[j]+=d[j]*center[j];}
   Eigen::LLT<MatrixXd> solver(G);if(solver.info()!=Eigen::Success){failed=true;return previous;}
   fit=solver.solve(rhs);if(solver.info()!=Eigen::Success){failed=true;return previous;}
 }
 for(int j=0;j<p;++j)previous[active[j]]=fit[j]/scale[active[j]];
 return previous;
}
VectorXd ridge(const Eigen::MappedSparseMatrix<double>& X,const VectorXd& z,const VectorXd& W,
 const std::vector<int>& active,const VectorXd& penalties,const VectorXd& prior,
 const VectorXd& scale,VectorXd previous,bool& failed,MatrixXd& cacheG,VectorXd& cacheXY,
 VectorXd& cacheW,bool cache,bool cache_weights) {
   auto rows=positive_rows(W);const int n=rows.size(),p=active.size();
  if(!n){failed=true;return previous;}
  VectorXd response(n),d(p),center(p),inv(p);
  for(int i=0;i<n;++i)response[i]=std::sqrt(W[rows[i]])*z[rows[i]];
  for(int j=0;j<p;++j) {
    int a=active[j];inv[j]=1/scale[a];d[j]=penalties[a]/(scale[a]*scale[a]);center[j]=prior[a]*scale[a];
  }
  SparseMatrix<double> A=weighted_design(X,rows,active,W);A=A*inv.asDiagonal();VectorXd fit(p);
if(p>n) {
   if(!L0glm::ridge_dual(A,response,d,center,fit,false)){failed=true;return previous;}
 }else {
   SparseMatrix<double> G=A.transpose()*A;VectorXd rhs=A.transpose()*response;
   for(int j=0;j<p;++j){G.coeffRef(j,j)+=d[j];rhs[j]+=d[j]*center[j];}
   G.makeCompressed();Eigen::SimplicialLLT<SparseMatrix<double>> solver;solver.compute(G);
   if(solver.info()!=Eigen::Success){failed=true;return previous;}
   fit=solver.solve(rhs);if(solver.info()!=Eigen::Success){failed=true;return previous;}
 }
 for(int j=0;j<p;++j)previous[active[j]]=fit[j]/scale[active[j]];
 return previous;
}
#include "glm_solver.h"
template<class Design>
void fit(const Design& X,State& s,const Family& f,const SolverControl& control) {
 const int n=X.rows(),p=X.cols();VectorXd old=s.beta;
 s.W=VectorXd::Zero(n);s.z=VectorXd::Zero(n);s.lambdas=VectorXd::Constant(p,R_PosInf);
  s.Wret=VectorXd::Zero(n);s.zret=VectorXd::Zero(n);
 for(int iteration=0;iteration<s.maxit;++iteration) {
  s.iterations=iteration+1;
  for(int i=0;i<n;++i) {
   double derivative=f.derivative(s.eta[i]),variance=f.variance(s.mu[i]),weight=0,response=0;
   if(f.scalar>=0) {
     const double scalar_mu=std::max(s.mu[i],L0glmScalar::tiny);
     auto a=L0glmScalar::derivatives(s.y[i],scalar_mu,f.scalar,f.parameter,f.theta,f.censor);
     double information=std::max({1e-20,a.information,s.y[i]==0?1/(f.parameter*scalar_mu):0.});
     weight=s.w[i]*information*derivative*derivative;
     response=s.eta[i]-s.offset[i]-a.gradient/(information*derivative);
   }else if(f.link==0){weight=s.w[i]/variance;response=s.y[i]-s.offset[i];}
   else if(f.kind==4&&f.link==1){weight=s.w[i]*std::exp(-s.eta[i]);response=s.eta[i]-s.offset[i]-1+s.y[i]/s.mu[i];}
   else {weight=s.w[i]*derivative*derivative/variance;response=s.eta[i]-s.offset[i]+(s.y[i]-s.mu[i])/derivative;}
   if(s.w[i]==0)weight=0;
    s.Wret[i]=weight;s.zret[i]=response;
   if(!(weight>0)||!std::isfinite(weight)||!std::isfinite(response)||derivative==0) {weight=0;response=0;}
   s.W[i]=weight;s.z[i]=response;
  }
  std::vector<int> active;
  for(int j=0;j<p;++j) {
   s.lambdas[j]=s.beta[j]!=0?s.lambda*s.phi/(4*s.beta[j]*s.beta[j]):R_PosInf;
   if(s.included[j])s.lambdas[j]=0;
   if(s.beta[j]!=0||s.included[j])active.push_back(j);
  }
  if(active.empty()){s.eta=s.offset;for(int i=0;i<n;++i)s.mu[i]=f.mean(s.eta[i]);s.converged=true;break;}
    const bool fast_llt=!control.external()&&control.method==0&&control.pivot==0&&control.ordering==0&&
      (!control.rowaug||std::is_same_v<Design,Eigen::MappedSparseMatrix<double>>);
    VectorXd candidate=fast_llt?ridge(X,s.z,s.W,active,s.lambdas,s.prior,s.scale,s.beta,s.solver_failed,
      s.cacheG,s.cacheXY,s.cacheW,s.cache,s.cache_weights):ridge_selected(X,s,active,control,iteration+1);
  for(int j=0;j<p;++j)candidate[j]=std::max(s.lower[j],std::min(s.upper[j],candidate[j]));
  bool support_changed=false;
  for(int j=0;j<p;++j)if(!s.included[j]&&((candidate[j]!=0)!=(old[j]!=0)))support_changed=true;
  VectorXd eta=predict(X,candidate,s.offset);
  if(s.free_intercept>=0) {
   double information=0,numerator=0;const int col=s.free_intercept;
   for(int i=0;i<n;++i)if(s.W[i]>0){double x=X.coeff(i,col);information+=s.W[i]*x*x;
     numerator+=s.W[i]*x*(s.z[i]-(eta[i]-s.offset[i]));}
   if(information>0&&std::isfinite(information)) {
    double a=numerator/information;candidate[col]+=a;
    for(int i=0;i<n;++i)eta[i]+=X.coeff(i,col)*a;
   }
  }
  auto penalty=[&](const VectorXd& b) {
    long double value=0;for(int j=0;j<p;++j)if(s.lambdas[j]>0&&std::isfinite(s.lambdas[j]))
      value+=s.lambdas[j]*(b[j]-s.prior[j])*(b[j]-s.prior[j]);return (double)value;
  };
  double reference=iteration>0&&!support_changed?s.previous_dev+penalty(old):R_PosInf;
  double threshold=std::isfinite(reference)?reference+1e-10*std::max(1.,std::abs(reference)):R_PosInf;
  VectorXd anchor=iteration==0?s.anchor:old,trial=candidate,trial_eta=eta;bool accepted=false;double dev=R_PosInf;
  for(int half=0;half<=50;++half) {
   dev=f.deviance(s.y,trial_eta,s.w);double value=dev+penalty(trial);
   if(std::isfinite(value)&&value<=threshold){accepted=true;break;}
   double fraction=std::ldexp(1.,-half-1);trial=anchor+fraction*(candidate-anchor);
   trial_eta=s.previous_eta+fraction*(eta-s.previous_eta);
  }
   if(!accepted){
     if(!std::isfinite(s.previous_dev))throw std::runtime_error("Cannot find finite IRLS means; supply feasible starts or offsets.");
     s.beta=anchor;s.eta=s.previous_eta;s.failed=true;
    for(int i=0;i<n;++i){s.mu[i]=f.mean(s.eta[i]);if(f.scalar<0&&(f.kind==1||f.kind==3||f.kind==4))s.mu[i]=std::max(eps,s.mu[i]);}break;}
  if(s.hard) {
   double largest=0;for(int j=0;j<p;++j)if(!s.included[j])largest=std::max(largest,std::abs(trial[j])*s.scale[j]);
   VectorXd cut=trial;
   bool changed=false;for(int j=0;j<p;++j)if(!s.included[j]&&trial[j]!=0&&
      (std::abs(trial[j])*s.scale[j]<largest*s.thresh||std::abs(trial[j])*s.scale[j]<1e-20)){cut[j]=0;changed=true;}
   if(changed) {
     VectorXd cut_eta=trial_eta+predict(X,VectorXd(cut-trial),VectorXd::Zero(n));double value=f.deviance(s.y,cut_eta,s.w);
    if(std::isfinite(value)){trial=cut;trial_eta=cut_eta;dev=value;}
   }
  }
  s.beta=trial;s.eta=trial_eta;
  for(int i=0;i<n;++i){s.mu[i]=f.mean(s.eta[i]);if(f.scalar<0&&(f.kind==1||f.kind==3||f.kind==4))s.mu[i]=std::max(eps,s.mu[i]);}
  support_changed=false;int selected=0;
  for(int j=0;j<p;++j){if(s.beta[j]!=0||s.included[j])++selected;
    if(!s.included[j]&&((s.beta[j]!=0)!=(old[j]!=0)))support_changed=true;}
  s.previous_dev=dev;
  if(s.estimate&&s.hard&&support_changed){if(s.frozen)++s.restarts;s.frozen=false;s.stable=0;}
  bool estimable=selected+s.offset_df<.7*s.npositive;
  if(!s.frozen&&s.hard&&iteration>=2&&estimable) {
   long double value=0;for(int i=0;i<n;++i)if(s.w[i]>0)
      value+=s.w[i]*(s.y[i]-s.mu[i])*(s.y[i]-s.mu[i])/f.variance(s.mu[i]);
   double target=(double)value/(s.npositive-selected-s.offset_df);
   if(std::isfinite(target)&&target>0) {
    bool stable=!support_changed&&std::abs(target/s.phi-1)<.01;s.stable=stable?s.stable+1:0;
    s.phi=std::min(2*s.phi,std::max(s.phi/2,target));++s.updates;s.frozen=s.stable>=3;
   }
  }
  double relative=0;
  if(!s.normalized)relative=(s.beta-old).norm()/std::max(1e-12,s.beta.norm());
  else {
   double numerator=0,denominator=0;bool included=false;
   for(int j=0;j<p;++j)if(s.included[j])included=true;else {
     numerator+=std::pow((s.beta[j]-old[j])*s.scale[j],2);denominator+=std::pow(s.beta[j]*s.scale[j],2);}
   relative=std::sqrt(numerator)/std::max(1e-12,std::sqrt(denominator));
   if(included&&s.old_prediction.size()) {
    VectorXd pred=s.eta-s.offset;double total=s.w.sum(),mean=total>0?s.w.dot(pred)/total:0;
    double change=std::sqrt((s.w.array()*(pred-s.old_prediction).array().square()).sum());
    double size=std::sqrt((s.w.array()*(pred.array()-mean).square()).sum());
    relative=std::max(relative,change/std::max(1e-12,std::max(std::sqrt(total),size)));
   }
  }
  if(relative<s.tol&&iteration>0&&(s.frozen||!s.hard||s.lambda==0||(iteration>=2&&!estimable))){s.converged=true;break;}
  old=s.beta;s.previous_eta=s.eta;if(s.old_prediction.size())s.old_prediction=s.eta-s.offset;
 }
}
void prefill_shared_gram(const Eigen::Map<MatrixXd>& X,std::vector<State>& states,const Family& family) {
 if(states.empty()||family.scalar>=0||!((family.kind==3&&family.link==1)||
     (family.kind==1&&family.link==2)))return;
 auto& first=states[0];const int p=X.cols();std::vector<int> good;
 bool ready=true;for(const auto& s:states)ready=ready&&s.cacheG.rows()==p;
 if(ready)return;
 for(int i=0;i<X.rows();++i)if(first.w[i]>0)good.push_back(i);
 const int n=good.size();if(!first.cache_weights||p>n||p>2000||!p)return;
 for(const auto& s:states)if(!s.cache_weights||s.cache||
    !(s.w.array()==first.w.array()).all()||!(s.scale.array()==first.scale.array()).all())return;
 const double multiplier=family.kind==1?4:1;
 MatrixXd A(n,p);VectorXd root(n);
 for(int i=0;i<n;++i)root[i]=std::sqrt(multiplier*first.w[good[i]]);
 for(int j=0;j<p;++j)for(int i=0;i<n;++i)
   A(i,j)=X(good[i],j)*root[i]/first.scale[j];
 MatrixXd G=MatrixXd::Zero(p,p);double one=1,zero=0;
 F77_CALL(dsyrk)("L","T",&p,&n,&one,A.data(),&n,&zero,G.data(),&p FCONE FCONE);
 for(auto& s:states) {
  s.cacheG=G;s.cacheXY=VectorXd::Zero(p);s.cacheW=multiplier*s.w;
 }
}
void prefill_shared_gram(const Eigen::MappedSparseMatrix<double>&,
 std::vector<State>&,const Family&) {}
template<class Design> std::vector<State> load_states(const Design& X,List contexts) {
  std::vector<State> states(contexts.size());
  for(int k=0;k<contexts.size();++k) {
    List a=contexts[k];State& s=states[k];
    s.y=as<VectorXd>(a["y"]);s.w=as<VectorXd>(a["weights"]);s.offset=as<VectorXd>(a["offset"]);
    s.beta=as<VectorXd>(a["beta"]);s.eta=as<VectorXd>(a["eta"]);s.mu=as<VectorXd>(a["mu"]);
    s.prior=as<VectorXd>(a["prior.mean"]);s.scale=as<VectorXd>(a["scaling"]);
    s.lower=as<VectorXd>(a["lower"]);s.upper=as<VectorXd>(a["upper"]);
    s.phi=as<double>(a["dispersion"]);s.lambda=as<double>(a["lambda"]);s.selection_lambda=s.lambda;
    s.tol=as<double>(a["tol"]);s.thresh=as<double>(a["thresh"]);s.maxit=as<int>(a["maxit"]);
    s.offset_df=as<int>(a["offset_df"]);s.estimate=as<bool>(a["estimate_dispersion"]);
    s.frozen=!s.estimate;s.normalized=as<bool>(a["normalized"]);s.hard=as<bool>(a["hard_threshold"]);
    s.cache=as<bool>(a["cache_gaussian"]);s.cache_weights=as<bool>(a["cache_fixed_weights"]);
    IntegerVector included=a["always.include"],intercept=a["free_intercept"];
    s.included.assign(X.cols(),0);s.allowed.assign(X.cols(),1);
    for(int i:included)if(i>0&&i<=X.cols())s.included[i-1]=1;
    s.free_intercept=intercept.size()?intercept[0]-1:-1;
    s.npositive=0;for(double w:s.w)if(w>0&&std::isfinite(w))++s.npositive;
    s.auto_beta=as<bool>(a["auto_beta"]);s.auto_phi=as<bool>(a["auto_phi"]);
    s.gaussian_name=as<bool>(a["gaussian_name"]);s.binomial_name=as<bool>(a["binomial_name"]);
    s.signal_family=as<bool>(a["signal_family"]);s.reference_p=as<double>(a["reference_p"]);
    if(a.containsElementNamed("allowed")) {
      LogicalVector allowed=a["allowed"];for(int j=0;j<X.cols();++j)s.allowed[j]=allowed[j];
    }
    s.center_shift=VectorXd::Zero(X.cols());s.original_prior=s.prior;
    if(a.containsElementNamed("center_shift")) {
      s.center_shift=as<VectorXd>(a["center_shift"]);s.center_intercept=as<int>(a["center_intercept"]);
      s.original_prior=as<VectorXd>(a["original_prior"]);
    }
    s.identity_boundary=as<bool>(a["identity_boundary"]);
  }
  return states;
}
#include "glm_group.h"
}
// [[Rcpp::export(rng = false)]]
List glm_native_cpp(SEXP design,List contexts,int family,int link,double theta=1,
 int scalar=-1,double parameter=1,double censor=0,int nthreads=1,
  NumericVector probes=NumericVector(),int max_passes=1,bool relaxed=false,bool nonnegative=false,
   List solver_control=List::create(),SEXP solver_bridge=R_NilValue,SEXP executor=R_NilValue) {
 L0glmNative::Family spec{family,link,scalar,theta,parameter,censor,
  -R::qnorm5(L0glmNative::eps,0,1,1,0),-R::qcauchy(L0glmNative::eps,0,1,1,0)};
 auto stream=as<std::vector<double>>(probes);
  L0glmNative::SolverControl control(solver_control,solver_bridge);
 if(Rf_isS4(design)) {
  Eigen::MappedSparseMatrix<double> X=as<Eigen::MappedSparseMatrix<double>>(design);
    return L0glmNative::model(X,contexts,spec,stream,nthreads,max_passes,relaxed,nonnegative,control,executor);
 }
 NumericMatrix mat(design);Eigen::Map<MatrixXd> X(mat.begin(),mat.nrow(),mat.ncol());
   return L0glmNative::model(X,contexts,spec,stream,nthreads,max_passes,relaxed,nonnegative,control,executor);
}
// [[Rcpp::export(rng = false)]]
List glm_native_pass_cpp(SEXP design,List contexts,int family,int link,double theta=1,
  int scalar=-1,double parameter=1,double censor=0,NumericVector probes=NumericVector(),
  List solver_control=List::create(),SEXP solver_bridge=R_NilValue) {
  L0glmNative::Family spec{family,link,scalar,theta,parameter,censor,
    -R::qnorm5(L0glmNative::eps,0,1,1,0),-R::qcauchy(L0glmNative::eps,0,1,1,0)};
  auto stream=as<std::vector<double>>(probes);L0glmNative::SolverControl control(solver_control,solver_bridge);
  if(Rf_isS4(design)) {
    Eigen::MappedSparseMatrix<double> X=as<Eigen::MappedSparseMatrix<double>>(design);
    return L0glmNative::channel_pass(X,contexts,spec,stream,control);
  }
  NumericMatrix mat(design);Eigen::Map<MatrixXd> X(mat.begin(),mat.nrow(),mat.ncol());
  return L0glmNative::channel_pass(X,contexts,spec,stream,control);
}
// [[Rcpp::export(rng = false)]]
List signal_start_cpp(SEXP design,NumericMatrix Y,NumericMatrix offset,NumericMatrix weights,
 NumericMatrix prior,IntegerVector mandatory,int nthreads=1) {
 if(Rf_isS4(design)) {
  Eigen::MappedSparseMatrix<double> X=as<Eigen::MappedSparseMatrix<double>>(design);
  return L0glmNative::signal_matrix(X,Y,offset,weights,prior,mandatory,nthreads);
 }
 NumericMatrix mat(design);Eigen::Map<MatrixXd> X(mat.begin(),mat.nrow(),mat.ncol());
 return L0glmNative::signal_matrix(X,Y,offset,weights,prior,mandatory,nthreads);
}
