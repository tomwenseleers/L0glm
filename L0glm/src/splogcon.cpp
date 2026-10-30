#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include "nuisance.h"
#include "scalar_families.h"
#include <algorithm>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace L0glmScalar;

namespace {
struct Parameters {
  int family, full_n, first;NumericVector phi, theta, weights;double censor;
  Parameters(int f, NumericVector p, NumericVector t, NumericVector w, int N, int a, double c):
    family(f), full_n(N), first(a), phi(p), theta(t), weights(w), censor(c) {}
  double p(int j) const {return phi[phi.size()==1 ? 0 : j];}
  double t(int j) const {return theta[theta.size()==1 ? 0 : j];}
  double w(int i,int j) const {return weights[weights.size()==1 ? 0 : first+i+full_n*j];}
};

// Bracketed exact Poisson coordinate solve; constant density terms cancel.
double poisson_coordinate(const std::vector<double>& y,const std::vector<double>& r,
 const std::vector<double>& x,const std::vector<double>& w,const std::vector<double>& phi,
 int family,double start,int maxit,double tol) {
  const int n=x.size();double hs=0,total=0,g0=0;
  for(int i=0;i<n;++i)if(x[i]>0&&w[i]>0) {
    const double weight=w[i]/(family==2?phi[i]:1);
    hs+=weight*x[i];total+=weight*y[i];g0+=weight*x[i]*(1-y[i]/std::max(r[i],tiny));
  }
  if(hs<=0)return 0;
  double fit=0;
  if(g0<0&&total>0) {
    double lo=0,hi=total/hs;
    fit=start>0?std::max(1e-12,std::min(start,hi)):.1*hi;
    for(int it=0;it<std::max(25,maxit);++it) {
      double g=hs,h=0;
      for(int i=0;i<n;++i)if(x[i]>0&&w[i]>0) {
        const double weight=w[i]/(family==2?phi[i]:1),mu=std::max(r[i]+fit*x[i],tiny);
        g-=weight*y[i]*x[i]/mu;h+=weight*y[i]*x[i]*x[i]/(mu*mu);
      }
      if(std::abs(g)<=tol*hs)break;
      if(g>0)hi=fit;else lo=fit;
      double next=h>0?fit-g/h:.5*(lo+hi);
      if(!(next>lo&&next<hi)||!std::isfinite(next))next=.5*(lo+hi);
      if(std::abs(next-fit)<=tol*std::max(1.,fit)){fit=next;break;}fit=next;
    }
  }
  auto delta=[&](double b) {long double d=0;
    for(int i=0;i<n;++i)if(x[i]>0&&w[i]>0)d+=w[i]*change(y[i],r[i]+start*x[i],r[i]+b*x[i],family,phi[i],1,.5);
    return (double)d;};
  double best=start,value=0;const double zero=delta(0),candidate=delta(fit);
  if(zero<value){best=0;value=zero;}if(candidate<value)best=fit;
  return best;
}

double fit_coordinate(const std::vector<double>& y,const std::vector<double>& r,
  const std::vector<double>& x,const std::vector<double>& w,const std::vector<double>& phi,
  const std::vector<double>& theta,int family,double censor,double start,int iterations,double tol,
  double ridge=0) {
  if((family==1||family==2)&&ridge==0) return poisson_coordinate(y,r,x,w,phi,family,start,iterations,tol);
  const int n=x.size();double b=start, prior=start;
  auto delta=[&](double old,double next) {
    long double value=0;
    for(int i=0;i<n;++i) if(x[i]>0 && w[i]>0)
      value+=w[i]*change(y[i],r[i]+x[i]*old,r[i]+x[i]*next,family,phi[i],theta[i],censor);
    value+=ridge*((next-prior)*(next-prior)-(old-prior)*(old-prior));
    return (double)value;
  };
  for(int it=0;it<iterations;++it) {
    long double g=ridge*(b-prior),h=ridge;
    for(int i=0;i<n;++i) if(x[i]>0 && w[i]>0) {
      const double mu=std::max(tiny,r[i]+b*x[i]);
      Terms a=derivatives(y[i],mu,family,phi[i],theta[i],censor);
      g+=w[i]*x[i]*a.gradient;h+=w[i]*x[i]*x[i]*a.information;
    }
    if(!(h>0) || !std::isfinite((double)g) || !std::isfinite((double)h)) break;
    const double proposal=std::max(0.,b-(double)(g/h));
    if(std::abs(proposal-b)<=tol*std::max(1.,b)) break;
    bool accepted=false;double next=proposal;
    for(int half=0;half<24;++half) {
      const double d=delta(b,next);
      if(std::isfinite(d) && d<=1e-10) {accepted=true;break;}
      next=.5*(next+b);
    }
    if(!accepted) break;
    const double distance=std::abs(next-b);b=next;
    if(distance<=tol*std::max(1.,b)) break;
  }
  if(b>0 && delta(b,0)<=0) b=0;
  return delta(start,b)<=1e-8 ? b : start;
}

struct State {NumericMatrix y,mu;State(NumericMatrix Y,NumericMatrix M):y(Y),mu(clone(M)) {}};
}

// [[Rcpp::export]]
SEXP slcn_state(NumericMatrix y,NumericMatrix mu) {return XPtr<State>(new State(y,mu),true);}
// [[Rcpp::export]]
List slcn_local(SEXP handle,int left,int right) {
  XPtr<State> s(handle);const int first=left-1,n=right-left+1,m=s->y.ncol();
  NumericMatrix y(n,m),mu(n,m);
  for(int j=0;j<m;++j)for(int i=0;i<n;++i){y(i,j)=s->y(first+i,j);mu(i,j)=s->mu(first+i,j);}
  y.attr("row_start")=left;return List::create(_["Y"]=y,_["mu"]=mu);
}
// [[Rcpp::export]]
void slcn_commit(SEXP handle,int left,NumericMatrix mu) {
  XPtr<State> s(handle);for(int j=0;j<mu.ncol();++j)for(int i=0;i<mu.nrow();++i)s->mu(left-1+i,j)=mu(i,j);
}
// [[Rcpp::export]]
NumericMatrix slcn_mean(SEXP handle) {XPtr<State> s(handle);return clone(s->mu);}

// [[Rcpp::export]]
double slcn_delta(NumericMatrix y,NumericMatrix old,NumericMatrix next_mu,int family,
  NumericVector phi,NumericVector theta,NumericVector weights,int full_n,int first=0,double censor=.5) {
  Parameters p(family,phi,theta,weights,full_n,first,censor);long double value=0;
  for(int j=0;j<y.ncol();++j)for(int i=0;i<y.nrow();++i)if(p.w(i,j)>0)
    value+=p.w(i,j)*change(y(i,j),old(i,j),next_mu(i,j),family,p.p(j),p.t(j),censor);
  return (double)value;
}

// [[Rcpp::export]]
double slcn_loglik(NumericMatrix y,NumericMatrix mu,int family,NumericVector phi,
  NumericVector theta,NumericVector weights,int nthreads=1,double censor=.5) {
  Parameters p(family,phi,theta,weights,y.nrow(),0,censor);NumericVector totals(y.ncol());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads)
#endif
  for(int j=0;j<y.ncol();++j) {
    long double value=0;
    for(int i=0;i<y.nrow();++i)if(p.w(i,j)>0)value-=p.w(i,j)*terms(y(i,j),mu(i,j),family,p.p(j),p.t(j),censor).loss;
    totals[j]=(double)value;
  }
  long double value=0;for(double a:totals)value+=a;return (double)value;
}

// Sparse profiles and background share one conditional nonnegative solver.
// Each thread owns a channel, avoiding per-component OpenMP barriers.
// [[Rcpp::export]]
List slcn_spectral(NumericMatrix y,NumericMatrix mean,S4 design,NumericMatrix spectra,
  int family,NumericVector phi,NumericVector theta,NumericVector weights,int full_n,int first=0,
  int sweeps=4,int nthreads=1,int maxit=12,double tol=1e-8,double censor=.5,bool keep_zeros=false) {
  IntegerVector ptr=design.slot("p"),index=design.slot("i"),dims=design.slot("Dim");
  NumericVector values=design.slot("x");const int n=y.nrow(),m=y.ncol(),k=dims[1];
  if(dims[0]!=n || spectra.nrow()!=m || spectra.ncol()!=k) stop("Incompatible design/spectrum dimensions.");
  Parameters par(family,phi,theta,weights,full_n,first,censor);
  NumericMatrix mu=clone(mean),S=clone(spectra);NumericVector changes(m),maximum_change(m);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads)
#endif
  for(int c=0;c<m;++c) {
    long double total=0;double biggest=0;
    std::vector<double> yy,r,h,w,pp,tt;
    for(int pass=0;pass<sweeps;++pass) {
      double pass_change=0;
      for(int j=0;j<k;++j) {
        const int len=ptr[j+1]-ptr[j];if(!len)continue;
        const double old=S(c,j);if(keep_zeros && old==0)continue;
        yy.resize(len);r.resize(len);h.resize(len);w.resize(len);pp.assign(len,par.p(c));tt.assign(len,par.t(c));
        for(int a=0;a<len;++a) {const int z=ptr[j]+a,i=index[z];h[a]=values[z];yy[a]=y(i,c);
          r[a]=std::max(0.,mu(i,c)-old*h[a]);w[a]=par.w(i,c);}
        const double next=fit_coordinate(yy,r,h,w,pp,tt,family,censor,old,maxit,tol);
        if(next==old)continue;
        long double d=0;
        for(int a=0;a<len;++a) {const int i=index[ptr[j]+a];const double candidate=r[a]+next*h[a];
          if(w[a]>0)d+=w[a]*change(yy[a],mu(i,c),candidate,family,pp[a],tt[a],censor);mu(i,c)=candidate;}
        total+=d;pass_change+=(double)d;S(c,j)=next;biggest=std::max(biggest,std::abs(next-old)/std::max(1.,old));
      }
      if(std::abs(pass_change)<1e-8)break;
    }
    changes[c]=(double)total;maximum_change[c]=biggest;
  }
  long double total=0;for(double a:changes)total+=a;
  return List::create(_["mu"]=mu,_["S"]=S,_["change"]=(double)total,
    _["channel_change"]=changes,_["maximum_change"]=maximum_change);
}

// Independent conditional row coefficients produce a free profile proposal;
// the R layer projects it to hard log concavity and checks the actual loss.
// [[Rcpp::export]]
NumericVector slcn_profile(NumericMatrix y,NumericMatrix mean,NumericVector old,
  NumericVector spectrum,int family,NumericVector phi,NumericVector theta,NumericVector weights,
  int full_n,int first=0,int maxit=8,int nthreads=1,double censor=.5,double ridge=0) {
  Parameters par(family,phi,theta,weights,full_n,first,censor);NumericVector out=clone(old);
  const int n=y.nrow(),m=y.ncol();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads)
#endif
  for(int i=0;i<n;++i) {
    std::vector<double> yy(m),r(m),x(m),w(m),pp(m),tt(m);
    for(int j=0;j<m;++j){yy[j]=y(i,j);x[j]=spectrum[j];r[j]=std::max(0.,mean(i,j)-old[i]*x[j]);
      w[j]=par.w(i,j);pp[j]=par.p(j);tt[j]=par.t(j);}
    out[i]=fit_coordinate(yy,r,x,w,pp,tt,family,censor,old[i],maxit,1e-9,ridge);
  }
  return out;
}

// One pooled log-theta score/curvature; all channels share the same theta.
// [[Rcpp::export]]
NumericVector slcn_nb_score(NumericMatrix y,NumericMatrix mu,double theta,
                           NumericVector weights,int nthreads=1) {
  const int n=y.nrow(),m=y.ncol();NumericMatrix result(m,2);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads)
#endif
  for(int j=0;j<m;++j) {
    long double U=0,H=0;
    for(int i=0;i<n;++i) {
      const double w=weights[weights.size()==1?0:i+n*j],Y=y(i,j),M=mu(i,j);
      if(w<=0 || (Y==0 && M==0))continue;
      U+=w*(R::digamma(Y+theta)-R::digamma(theta)+std::log(theta)+1-std::log(theta+M)-(Y+theta)/(theta+M));
      H+=w*(R::trigamma(Y+theta)-R::trigamma(theta)+1/theta-2/(theta+M)+(Y+theta)/((theta+M)*(theta+M)));
    }
    result(j,0)=(double)U;result(j,1)=(double)H;
  }
  long double U=0,H=0;for(int j=0;j<m;++j){U+=result(j,0);H+=result(j,1);}
  return NumericVector::create(_["score"]=theta*(double)U,
    _["curvature"]=theta*theta*(double)H+theta*(double)U);
}

// Joint constrained IRLS for the small spline block. Sparse row bookkeeping
// avoids treating the local B-spline basis as a dense n by q operation.
// [[Rcpp::export]]
List slcn_background(NumericMatrix y,NumericMatrix signal,NumericMatrix basis,NumericMatrix start,
  int family,NumericVector phi,NumericVector theta,NumericVector weights,int nthreads=1,
  int maxit=40,double tol=1e-8,double censor=.5) {
  using Eigen::MatrixXd;using Eigen::VectorXd;
  const int n=y.nrow(),m=y.ncol(),q=basis.ncol();
  if(basis.nrow()!=n || signal.nrow()!=n || signal.ncol()!=m || start.nrow()!=q || start.ncol()!=m)
    stop("Incompatible background dimensions.");
  Eigen::Map<const MatrixXd> B(basis.begin(),n,q);
  NumericMatrix coefficients=clone(start),background(n,m);NumericVector iterations(m),changes(m);
  Parameters par(family,phi,theta,weights,n,0,censor);
  std::vector<int> ptr(n+1),index;std::vector<double> value;
  for(int i=0;i<n;++i){for(int j=0;j<q;++j)if(B(i,j)!=0){index.push_back(j);value.push_back(B(i,j));}ptr[i+1]=index.size();}
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads)
#endif
  for(int channel=0;channel<m;++channel) {
    VectorXd a(q),fixed(n);for(int j=0;j<q;++j)a[j]=start(j,channel);
    for(int i=0;i<n;++i)fixed[i]=signal(i,channel);
    VectorXd mu=fixed+B*a;long double total=0;
    auto delta=[&](const VectorXd& next){long double d=0;
      for(int i=0;i<n;++i)if(par.w(i,channel)>0)d+=par.w(i,channel)*
        change(y(i,channel),mu[i],next[i],family,par.p(channel),par.t(channel),censor);
      return (double)d;};
    for(int iteration=0;iteration<maxit;++iteration) {
      MatrixXd G=MatrixXd::Zero(q,q);VectorXd gradient=VectorXd::Zero(q);
      for(int i=0;i<n;++i)if(par.w(i,channel)>0) {
        const Terms d=derivatives(y(i,channel),std::max(tiny,mu[i]),family,par.p(channel),par.t(channel),censor);
        const double g=par.w(i,channel)*d.gradient,h=par.w(i,channel)*std::max(0.,d.information);
        for(int j=ptr[i];j<ptr[i+1];++j) {
          const int u=index[j];gradient[u]+=value[j]*g;
          for(int k=ptr[i];k<ptr[i+1];++k)G(u,index[k])+=value[j]*value[k]*h;
        }
      }
      if(!G.allFinite() || !gradient.allFinite())break;
      const VectorXd proposal=L0glm::nuisance_qp(G,G*a-gradient,std::vector<char>(q,1));
      if(!proposal.allFinite())break;
      VectorXd next=a,pred=mu;bool accepted=false;double fraction=1,change_value=0;
      for(int half=0;half<25;++half) {
        next=a+fraction*(proposal-a);pred=fixed+B*next;change_value=delta(pred);
        if(std::isfinite(change_value) && change_value<=1e-8){accepted=true;break;}
        fraction*=.5;
      }
      if(!accepted)break;
      const double distance=(next-a).norm()/std::max(1.,a.norm());
      a=next;mu=pred;total+=change_value;iterations[channel]=iteration+1;
      for(int j=0;j<q;++j)if(a[j]>0) {
        pred=mu-B.col(j)*a[j];for(int i=0;i<n;++i)pred[i]=std::max(0.,pred[i]);
        const double d=delta(pred);
        if(std::isfinite(d) && d<=1e-10){a[j]=0;mu=pred;total+=d;}
      }
      double kkt=0;
      for(int j=0;j<q;++j)kkt=std::max(kkt,std::abs(a[j]>0 ? gradient[j] : std::min(0.,gradient[j])));
      kkt/=std::max(1.,(G*a).lpNorm<Eigen::Infinity>());
      if(distance<tol && kkt<std::sqrt(tol))break;
    }
    VectorXd fitted=B*a;
    for(int j=0;j<q;++j)coefficients(j,channel)=a[j];
    for(int i=0;i<n;++i)background(i,channel)=fitted[i];
    changes[channel]=(double)total;
  }
  return List::create(_["coefficients"]=coefficients,_["background"]=background,
    _["iterations"]=iterations,_["channel_change"]=changes);
}
