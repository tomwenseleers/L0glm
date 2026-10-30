#include "scalar_families.h"
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace L0glmScalar;

namespace {
double mean_mle(const std::vector<double>& y,const std::vector<double>& w,
                int f,double phi,double censor) {
  double total=0,mass=0,maximum=0;
  for(size_t i=0;i<y.size();++i) if(w[i]>0) {
    total+=w[i]*y[i];mass+=w[i];maximum=std::max(maximum,y[i]);
  }
  if(!mass || !maximum) return 0;
  if(f==5 && phi==1) return total/mass;
  if(f==4 && std::all_of(y.begin(),y.end(),[](double x){return x>0;})) return total/mass;
  auto score=[&](double mu) {
    double g=0;
    for(size_t i=0;i<y.size();++i) if(w[i]>0)
      g+=w[i]*derivatives(y[i],mu,f,phi,1,censor).gradient;
    return g;
  };
  double lo=0,hi=std::max({1.,maximum,phi}),b=std::max(1e-12,total/mass);
  for(int it=0;it<100 && score(hi)<0;++it) hi*=2;
  b=std::min(b,hi);
  for(int it=0;it<80;++it) {
    double g=0,h=0;
    for(size_t i=0;i<y.size();++i) if(w[i]>0) {
      const Terms d=derivatives(y[i],std::max(b,1e-300),f,phi,1,censor);
      g+=w[i]*d.gradient;h+=w[i]*std::max(0.,d.information);
    }
    if(g>0) hi=b;else lo=b;
    double next=h>0 ? b-g/h : .5*(lo+hi);
    if(!std::isfinite(next) || next<=lo || next>=hi) next=.5*(lo+hi);
    if(hi-lo<=1e-11*std::max(1.,b)) return .5*(lo+hi);
    if(std::abs(g)<=1e-11*mass/std::max(1.,b)) return b;
    b=next;
  }
  return b;
}
double saturated(double y,int f,double phi,double censor) {
  if(y==0) return 0;
  if(f==4 || (f==5 && phi==1)) return y;
  return mean_mle({y},{1},f,phi,censor);
}
double weight_at(const NumericVector& w,int i,int j,int n) {
  return w[w.size()==1?0:w.size()==n?i:i+n*j];
}
void parameters(int f,double phi,double censor) {
  if(f<4||f>6||!std::isfinite(phi)||phi<=0||!std::isfinite(censor)||censor<0||
      (f==5&&phi<1)) stop("Invalid scalar family parameters.");
}
void response(double y,int f,double censor) {
  if(!std::isfinite(y)||y<0||(f==5&&std::abs(y-std::round(y))>1e-7)||
      (f!=5&&(y==0?censor==0:y<censor))) stop("Response outside the scalar family domain.");
}
}

// [[Rcpp::export]]
NumericMatrix L0glm_scalar_terms_cpp(NumericVector y,NumericVector mu,int family,
                                    double phi,double censor,bool deviance=false) {
  parameters(family,phi,censor);
  if(y.size()!=mu.size()) stop("Response and mean lengths differ.");
  NumericMatrix out(y.size(),deviance?4:3);
  for(R_xlen_t i=0;i<y.size();++i) {
    response(y[i],family,censor);
    const Terms t=terms(y[i],mu[i],family,phi,1,censor);
    const Terms d=derivatives(y[i],std::max(mu[i],tiny),family,phi,1,censor);
    out(i,0)=-t.loss;out(i,1)=d.gradient;
    out(i,2)=std::max({1e-20,d.information,
      y[i]==0?1/(phi*std::max(mu[i],tiny)):0.});
    if(deviance) out(i,3)=std::max(0.,2*(t.loss-terms(y[i],
      saturated(y[i],family,phi,censor),family,phi,1,censor).loss));
  }
  return out;
}

// Exact likelihood-ratio pooling for scalar families whose MLE is not y.
// [[Rcpp::export]]
List L0glm_scalar_profile_cpp(NumericMatrix y,NumericVector weights,IntegerMatrix labels,
                              int family,double phi,double censor,int nthreads=1) {
  parameters(family,phi,censor);
  const int n=y.nrow(),m=y.ncol(),cuts=labels.ncol();
  if(labels.nrow()!=m||!(weights.size()==1||weights.size()==n||weights.size()==n*m)||nthreads<1)
    stop("Invalid profile dimensions or thread count.");
  for(double v:y) response(v,family,censor);
  for(double v:weights) if(!std::isfinite(v)||v<0) stop("Invalid weights.");
  for(int g:labels) if(g<1||g>m) stop("Invalid group label.");
  NumericVector dev(cuts),df(cuts);NumericMatrix centres(n,m);
  for(int cut=0;cut<cuts;++cut) {
    int k=0;for(int j=0;j<m;++j)k=std::max(k,labels(j,cut));
    NumericVector row_loss(n),row_df(n);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(nthreads)
    #endif
    for(int i=0;i<n;++i) {
      std::vector<std::vector<double>> yy(k),ww(k);
      for(int j=0;j<m;++j) {int g=labels(j,cut)-1;yy[g].push_back(y(i,j));ww[g].push_back(weight_at(weights,i,j,n));}
      double loss=0,count=0;
      for(int g=0;g<k;++g) {
        const double mu=mean_mle(yy[g],ww[g],family,phi,censor);
        bool used=false;
        for(size_t a=0;a<yy[g].size();++a) if(ww[g][a]>0) {
          used=true;loss+=2*ww[g][a]*(terms(yy[g][a],mu,family,phi,1,censor).loss-
            terms(yy[g][a],saturated(yy[g][a],family,phi,censor),family,phi,1,censor).loss);
        }
        if(used)++count;
        if(cut==cuts-1) centres(i,g)=mu;
      }
      row_loss[i]=loss;row_df[i]=count;
    }
    long double sum=0;for(int i=0;i<n;++i){sum+=row_loss[i];df[cut]+=row_df[i];}
    dev[cut]=std::max(0.,(double)sum);
  }
  int k=0;for(int j=0;j<m;++j)k=std::max(k,labels(j,cuts-1));
  return List::create(_["deviance"]=dev,_["df"]=df,_["centres"]=centres(_,Range(0,k-1)));
}

// [[Rcpp::export]]
NumericVector L0glm_scalar_dist_cpp(NumericMatrix y,NumericVector weights,int family,
                                   double phi,double censor,int nthreads=1) {
  parameters(family,phi,censor);const int n=y.nrow(),m=y.ncol();
  if(!(weights.size()==1||weights.size()==n||weights.size()==n*m)||nthreads<1) stop("Invalid distance controls.");
  for(double v:y)response(v,family,censor);
  for(double v:weights)if(!std::isfinite(v)||v<0)stop("Invalid weights.");
  NumericVector out((R_xlen_t)m*(m-1)/2);std::vector<int> left,right;
  for(int a=0;a<m;++a)for(int b=a+1;b<m;++b){left.push_back(a);right.push_back(b);}
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(nthreads)
  #endif
  for(R_xlen_t pair=0;pair<out.size();++pair) {
    const int a=left[pair],b=right[pair];long double loss=0;
    for(int i=0;i<n;++i) {
      std::vector<double> yy={y(i,a),y(i,b)},ww={weight_at(weights,i,a,n),weight_at(weights,i,b,n)};
      const double mu=mean_mle(yy,ww,family,phi,censor);
      for(int j=0;j<2;++j)if(ww[j]>0) loss+=2*ww[j]*(terms(yy[j],mu,family,phi,1,censor).loss-
        terms(yy[j],saturated(yy[j],family,phi,censor),family,phi,1,censor).loss);
    }
    out[pair]=std::max(0.,(double)loss);
  }
  return out;
}
