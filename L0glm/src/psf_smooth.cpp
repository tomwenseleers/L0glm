#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "scalar_families.h"
using namespace Rcpp;

namespace {
std::vector<double> concave_majorant(const NumericVector& value) {
  const int n=value.size();std::vector<int> knots;std::vector<double> out(n);
  for(int i=0;i<n;++i) {
    knots.push_back(i);
    while(knots.size()>2) {
      const int k=knots.size(),a=knots[k-3],b=knots[k-2],c=knots[k-1];
      if((value[b]-value[a])/(b-a)>=(value[c]-value[b])/(c-b))break;
      knots[k-2]=c;knots.pop_back();
    }
  }
  for(size_t j=1;j<knots.size();++j) {
    const int a=knots[j-1],b=knots[j];
    for(int i=a;i<b;++i)out[i]=value[a]+(value[b]-value[a])*(i-a)/(b-a);
  }
  if(n)out[n-1]=value[n-1];return out;
}
}

// Weighted log-profile least squares with roughness and concavity constraints.
// ADMM uses one banded Cholesky factorization; every iteration is linear in n.
// [[Rcpp::export]]
List logconc_quadratic_cpp(NumericVector target,NumericVector weights,double smooth=0,
                    int maxit=2000,double tol=1e-8) {
  const int n=target.size(),q=n-2;
  if(n<3||weights.size()!=n||!std::isfinite(smooth)||smooth<0||maxit<1||tol<=0)
    stop("Invalid log-concave smoothing controls.");
  double wsum=0;
  for(int i=0;i<n;++i) {
    if(!std::isfinite(target[i])||!std::isfinite(weights[i])||weights[i]<=0)
      stop("Finite working targets and positive weights are required.");
    wsum+=weights[i];
  }
  std::vector<double> w(n),v(n),rhs(n),z(n),last(n),u(q),dual(q,0),uold(q),
    diag(n),sub1(n,0),sub2(n,0),ld(n),l1(n),l2(n),forward(n),box(n),box_dual(n,0);
  const double rho=1,penalty=smooth+rho;
  for(int i=0;i<n;++i) {
    w[i]=weights[i]*n/wsum;v[i]=target[i];
    z[i]=v[i];diag[i]=w[i]+rho;box[i]=std::max(-600.,v[i]);
  }
  for(int j=0;j<q;++j) {
    diag[j]+=penalty;diag[j+1]+=4*penalty;diag[j+2]+=penalty;
    sub1[j+1]-=2*penalty;sub1[j+2]-=2*penalty;sub2[j+2]+=penalty;
    u[j]=std::min(0.,v[j]-2*v[j+1]+v[j+2]);
  }
  for(int i=0;i<n;++i) {
    l2[i]=i>=2?sub2[i]/ld[i-2]:0;
    l1[i]=i>=1?(sub1[i]-(i>=2?l2[i]*l1[i-1]:0))/ld[i-1]:0;
    ld[i]=std::sqrt(diag[i]-l1[i]*l1[i]-l2[i]*l2[i]);
    if(!std::isfinite(ld[i])||ld[i]<=0)stop("Smoothing factorization failed.");
  }
  int iteration=0;double primal=0,residual=0;
  for(;iteration<maxit;++iteration) {
    last=z;uold=u;
    for(int i=0;i<n;++i)rhs[i]=w[i]*v[i]+rho*(box[i]-box_dual[i]);
    for(int j=0;j<q;++j) {
      const double a=rho*(u[j]-dual[j]);rhs[j]+=a;rhs[j+1]-=2*a;rhs[j+2]+=a;
    }
    for(int i=0;i<n;++i)forward[i]=(rhs[i]-(i?l1[i]*forward[i-1]:0)-
      (i>=2?l2[i]*forward[i-2]:0))/ld[i];
    for(int i=n-1;i>=0;--i)z[i]=(forward[i]-(i+1<n?l1[i+1]*z[i+1]:0)-
      (i+2<n?l2[i+2]*z[i+2]:0))/ld[i];
    primal=0;residual=0;
    for(int i=0;i<n;++i) {
      const double old=box[i];box[i]=std::max(-600.,z[i]+box_dual[i]);
      box_dual[i]+=z[i]-box[i];primal=std::max(primal,std::abs(z[i]-box[i]));
      residual=std::max(residual,rho*std::abs(box[i]-old));
    }
    for(int j=0;j<q;++j) {
      const double d=z[j]-2*z[j+1]+z[j+2];u[j]=std::min(0.,d+dual[j]);
      dual[j]+=d-u[j];primal=std::max(primal,std::abs(d-u[j]));
      residual=std::max(residual,rho*std::abs(u[j]-uold[j]));
    }
    if(primal<tol&&residual<tol){++iteration;break;}
    if(iteration%100==0)checkUserInterrupt();
  }
  NumericVector h(n),logfit(n);double top=*std::max_element(z.begin(),z.end()),sum=0;
  for(int i=0;i<n;++i){h[i]=std::exp(z[i]-top);logfit[i]=z[i];sum+=h[i];}
  for(int i=0;i<n;++i)h[i]/=sum;
  return List::create(_["h"]=h,_["logfit"]=logfit,_["iterations"]=iteration,
    _["converged"]=(primal<tol&&residual<tol),_["primal"]=primal,_["dual"]=residual);
}

// Original observation likelihood; the exponential parameterizes the profile.
// [[Rcpp::export]]
List logconc_profile_cpp(NumericMatrix y,NumericMatrix mean,NumericVector profile,
  NumericVector spectrum,int family,NumericVector phi,NumericVector theta,
  NumericVector weights,int full_n,int first=0,int maxit=12,double censor=.5,
  NumericVector ridge=NumericVector::create(0),NumericVector prior=NumericVector::create(0)) {
  const int n=y.nrow(),m=y.ncol();
  if(mean.nrow()!=n||mean.ncol()!=m||profile.size()!=n||spectrum.size()!=m||
     full_n<first+n||first<0||maxit<1||
     !(weights.size()==1||weights.size()==full_n*m)||
     !(phi.size()==1||phi.size()==m)||!(theta.size()==1||theta.size()==m)||
     !(ridge.size()==1||ridge.size()==n)||!(prior.size()==1||prior.size()==n))
    stop("Invalid log-concave profile dimensions.");
  auto weight=[&](int i,int j){return weights[weights.size()==1?0:first+i+full_n*j];};
  auto dispersion=[&](int j){return phi[phi.size()==1?0:j];};
  auto shape=[&](int j){return theta[theta.size()==1?0:j];};
  NumericVector h=clone(profile),eta(n),information(n),target(n);
  for(int i=0;i<n;++i)if(!std::isfinite(h[i])||h[i]<=0)stop("A profile needs positive finite support.");
  for(int j=0;j<m;++j)if(!std::isfinite(spectrum[j])||spectrum[j]<0)stop("Spectra must be nonnegative.");
  auto objective=[&](const NumericVector& value) {
    long double total=0;
    for(int i=0;i<n;++i) {
      for(int j=0;j<m;++j)if(weight(i,j)>0) {
        const double mu=std::max(0.,mean(i,j)-profile[i]*spectrum[j])+value[i]*spectrum[j];
        total+=weight(i,j)*L0glmScalar::terms(y(i,j),mu,family,dispersion(j),shape(j),censor).loss;
      }
      const double a=ridge[ridge.size()==1?0:i],b=prior[prior.size()==1?0:i];
      total+=.5*a*(value[i]-b)*(value[i]-b);
    }
    return (double)total;
  };
  const double initial=objective(h);double loss=initial;int iteration=0;
  std::vector<double> history;
  for(;iteration<maxit;++iteration) {
    double maximum=0;
    NumericVector gradient(n);
    for(int i=0;i<n;++i) {
      eta[i]=std::log(h[i]);long double g=0,I=0;
      for(int j=0;j<m;++j)if(weight(i,j)>0&&spectrum[j]>0) {
        const double a=h[i]*spectrum[j];
        const double mu=std::max(0.,mean(i,j)-profile[i]*spectrum[j])+a;
        const auto d=L0glmScalar::derivatives(y(i,j),std::max(mu,L0glmScalar::tiny),
          family,dispersion(j),shape(j),censor);
        g+=weight(i,j)*a*d.gradient;I+=weight(i,j)*a*a*d.information;
      }
      const double a=ridge[ridge.size()==1?0:i],b=prior[prior.size()==1?0:i];
      gradient[i]=(double)g+a*h[i]*(h[i]-b);
      information[i]=(double)I+a*h[i]*h[i];maximum=std::max(maximum,information[i]);
    }
    if(!(maximum>0))break;
    for(int i=0;i<n;++i) {
      information[i]+=std::max(1e-12,maximum*1e-7);
      target[i]=eta[i]-gradient[i]/information[i];
    }
    NumericVector proposed;
    if(n>=3) {
      List fit=logconc_quadratic_cpp(target,information,0,2000,1e-7);
      proposed=fit["logfit"];
    } else proposed=target;
    const auto feasible=concave_majorant(proposed);
    NumericVector candidate(n);bool accepted=false;double change=0;
    for(int half=0;half<=25;++half) {
      const double fraction=std::ldexp(1.,-half);
      bool finite=true;
      for(int i=0;i<n;++i) {candidate[i]=std::exp(eta[i]+fraction*(feasible[i]-eta[i]));
        finite=finite&&std::isfinite(candidate[i])&&candidate[i]>0;}
      if(!finite)continue;
      const double value=objective(candidate);
      if(std::isfinite(value)&&value<=loss+1e-9) {
        for(int i=0;i<n;++i)change=std::max(change,std::abs(candidate[i]-h[i]));
        h=candidate;loss=value;accepted=true;break;
      }
    }
    if(!accepted)break;
    history.push_back(loss);double size=0;for(double a:h)size=std::max(size,a);
    if(change<1e-8*std::max(1.,size)){++iteration;break;}
    checkUserInterrupt();
  }
  return List::create(_["profile"]=h,_["initial_loss"]=initial,_["loss"]=loss,
    _["iterations"]=iteration,_["history"]=wrap(history));
}
