#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

namespace {
struct SpoloState {
  NumericMatrix y, mean;
  SpoloState(NumericMatrix observed, NumericMatrix mu): y(observed), mean(clone(mu)) {}
};
void finite_nonnegative(const NumericMatrix& x, const char* name) {
  for(R_xlen_t i=0;i<x.length();++i)
    if(!std::isfinite(x[i]) || x[i]<0) stop("%s must be finite and nonnegative.",name);
}
void same_dimensions(const NumericMatrix& a,const NumericMatrix& b) {
  if(a.nrow()!=b.nrow() || a.ncol()!=b.ncol()) stop("Incompatible Poisson matrix dimensions.");
}
double poisson_change(double y,double old_mu,double next_mu) {
  if(!std::isfinite(next_mu)||next_mu<0)return R_PosInf;
  if(next_mu==old_mu)return 0;
  if(y==0)return 2*(next_mu-old_mu);
  if(next_mu<=0)return R_PosInf;
  if(old_mu<=0)return R_NegInf;
  const double change=next_mu-old_mu,ratio=change/old_mu;
  const double log_ratio=std::abs(ratio)<.5?std::log1p(ratio):std::log(next_mu)-std::log(old_mu);
  return 2*(change-y*log_ratio);
}
double coordinate_mle(const double* y,const double* r,const double* h,int n,
                      double start,int maxit,double tol) {
  double hs=0,total=0,g0=0;
  for(int i=0;i<n;++i)if(h[i]>0){hs+=h[i];total+=y[i];g0+=h[i]*(1-y[i]/std::max(r[i],1e-15));}
  if(hs<=0)return 0;
  double fit=0;
  if(g0<0&&total>0) {
    double lo=0,hi=total/hs;
    fit=std::max(1e-12,std::min(start,hi));if(start<=0)fit=.1*hi;
    for(int it=0;it<maxit;++it) {
      double g=hs,w=0;
      for(int i=0;i<n;++i)if(h[i]>0) {
        const double mu=std::max(r[i]+fit*h[i],1e-15);
        g-=y[i]*h[i]/mu;w+=y[i]*h[i]*h[i]/(mu*mu);
      }
      if(std::abs(g)<=tol*hs)break;
      if(g>0)hi=fit;else lo=fit;
      double next=w>0?fit-g/w:.5*(lo+hi);
      if(!(next>lo&&next<hi)||!std::isfinite(next))next=.5*(lo+hi);
      if(std::abs(next-fit)<=tol*(1+fit)){fit=next;break;}fit=next;
    }
  }
  auto change=[&](double b) {
    long double value=0;
    for(int i=0;i<n;++i)if(h[i]>0)value+=poisson_change(y[i],r[i]+start*h[i],r[i]+b*h[i]);
    return (double)value;
  };
  double best=start,value=0;const double zero=change(0);
  if(zero<value){best=0;value=zero;}if(change(fit)<value)best=fit;
  return best;
}
double js_distance(const double* spectrum,double area,const double* p,int m) {
  double value=0;
  for(int c=0;c<m;++c) {
    const double q=spectrum[c]/area,mid=.5*(q+p[c]);
    if(q>0)value+=.5*q*std::log(q/mid);if(p[c]>0)value+=.5*p[c]*std::log(p[c]/mid);
  }
  return value;
}
}

// [[Rcpp::export]]
SEXP spolo_state_cpp(NumericMatrix y,NumericMatrix mu) {
  same_dimensions(y,mu);finite_nonnegative(y,"y");finite_nonnegative(mu,"mu");
  return XPtr<SpoloState>(new SpoloState(y,mu),true);
}
// [[Rcpp::export]]
List spolo_local_cpp(SEXP handle,int left,int right) {
  XPtr<SpoloState> s(handle);
  if(left<1||right<left||right>s->y.nrow())stop("Invalid local profile interval.");
  --left;const int n=right-left,m=s->y.ncol();NumericMatrix y(n,m),mu(n,m);
  for(int c=0;c<m;++c)for(int i=0;i<n;++i){y(i,c)=s->y(left+i,c);mu(i,c)=s->mean(left+i,c);}
  return List::create(_["Y"]=y,_["mu"]=mu);
}
// [[Rcpp::export]]
void spolo_commit_cpp(SEXP handle,int left,NumericMatrix mu) {
  XPtr<SpoloState> s(handle);
  if(left<1||left+mu.nrow()-1>s->mean.nrow()||mu.ncol()!=s->mean.ncol())stop("Invalid local mean update.");
  --left;for(int c=0;c<mu.ncol();++c)for(int i=0;i<mu.nrow();++i)s->mean(left+i,c)=mu(i,c);
}
// [[Rcpp::export]]
NumericMatrix spolo_mean_cpp(SEXP handle) {XPtr<SpoloState> s(handle);return clone(s->mean);}
// [[Rcpp::export]]
double spolo_delta_cpp(NumericMatrix y,NumericMatrix old_mu,NumericMatrix new_mu) {
  same_dimensions(y,old_mu);same_dimensions(y,new_mu);long double value=0;
  for(R_xlen_t i=0;i<y.length();++i)value+=poisson_change(y[i],old_mu[i],new_mu[i]);
  return (double)value;
}

// Conditional nonnegative Poisson MLEs. Sparse profile columns are read directly;
// full chromatograms never require a dense n-scans by n-compounds design.
// [[Rcpp::export]]
List spolo_spectral_cpp(NumericMatrix y,NumericMatrix mean,SEXP profiles,
                        NumericMatrix spectra,int sweeps=1,int nthreads=1,
                        int maxit=25,double tol=1e-7) {
  same_dimensions(y,mean);finite_nonnegative(y,"y");finite_nonnegative(mean,"mean");
  finite_nonnegative(spectra,"spectra");
  if(sweeps<0||nthreads<1||maxit<1||!std::isfinite(tol)||tol<=0)stop("Invalid spectral solver controls.");
  const int n=y.nrow(),m=y.ncol(),k=spectra.ncol();
  bool sparse=Rf_isS4(profiles)&&Rf_inherits(profiles,"dgCMatrix");
  NumericMatrix dense;IntegerVector ptr,index;NumericVector values;
  if(sparse) {
    S4 H(profiles);IntegerVector dims=H.slot("Dim");
    if(dims[0]!=n||dims[1]!=k)stop("Profiles and spectra have incompatible dimensions.");
    ptr=H.slot("p");index=H.slot("i");values=H.slot("x");
    for(double v:values)if(!std::isfinite(v)||v<0)stop("Profiles must be finite and nonnegative.");
  } else {
    if(Rf_isS4(profiles)){Function as_matrix=Environment::base_env()["as.matrix"];dense=as<NumericMatrix>(as_matrix(profiles));}
    else dense=as<NumericMatrix>(profiles);
    if(dense.nrow()!=n||dense.ncol()!=k)stop("Profiles and spectra have incompatible dimensions.");
    finite_nonnegative(dense,"profiles");
  }
  if(spectra.nrow()!=m)stop("Spectra must have one row per outcome channel.");
  NumericMatrix mu=clone(mean),S=clone(spectra);NumericVector history(sweeps);double total_change=0;
  for(int sweep=0;sweep<sweeps;++sweep) {
    double sweep_change=0;
    for(int j=0;j<k;++j) {
      if(j%32==0)checkUserInterrupt();
      int left=n,right=-1;
      if(sparse){for(int z=ptr[j];z<ptr[j+1];++z)if(values[z]>0){left=std::min(left,index[z]);right=std::max(right,index[z]);}}
      else {left=0;right=n-1;while(left<n&&dense(left,j)<=0)++left;while(right>=left&&dense(right,j)<=0)--right;}
      if(right<left)continue;
      const int len=right-left+1;std::vector<double> h(len,0);
      if(sparse){for(int z=ptr[j];z<ptr[j+1];++z)if(index[z]>=left&&index[z]<=right)h[index[z]-left]=values[z];}
      else for(int i=0;i<len;++i)h[i]=dense(left+i,j);
      std::vector<double> changes(m,0);
      #ifdef _OPENMP
      #pragma omp parallel for schedule(static) num_threads(nthreads)
      #endif
      for(int c=0;c<m;++c) {
        std::vector<double> observed(len),r(len);const double old=S(c,j);
        for(int i=0;i<len;++i){observed[i]=y(left+i,c);r[i]=std::max(0.0,mu(left+i,c)-h[i]*old);}
        const double next=coordinate_mle(observed.data(),r.data(),h.data(),len,old,maxit,tol);
        long double change=0;
        for(int i=0;i<len;++i){const double value=r[i]+h[i]*next;change+=poisson_change(observed[i],mu(left+i,c),value);mu(left+i,c)=value;}
        S(c,j)=next;changes[c]=(double)change;
      }
      for(double v:changes)sweep_change+=v;
    }
    history[sweep]=sweep_change;total_change+=sweep_change;
  }
  return List::create(_["mu"]=mu,_["S"]=S,_["change"]=total_change,_["history"]=history);
}

// [[Rcpp::export]]
IntegerVector spolo_connected_cpp(NumericMatrix spectra,NumericVector area,
  NumericVector center,IntegerVector lo,IntegerVector hi,int anchor,double radius) {
  const int m=spectra.nrow(),k=spectra.ncol();--anchor;
  if(area.size()!=k||lo.size()!=k||hi.size()!=k||center.size()!=m)stop("Invalid spectral-clustering dimensions.");
  if(anchor<0||anchor>=k||area[anchor]<=0)return IntegerVector();
  std::vector<int> state(k,0);state[anchor]=2;int left=lo[anchor],right=hi[anchor];bool changed=true;
  while(changed) {
    checkUserInterrupt();changed=false;
    for(int j=0;j<k;++j) {
      if(state[j]||area[j]<=0||lo[j]>right||hi[j]<left)continue;
      const double distance=js_distance(&spectra(0,j),area[j],center.begin(),m);
      state[j]=distance<=radius?2:1;
      if(state[j]==2){left=std::min(left,lo[j]);right=std::max(right,hi[j]);changed=true;}
    }
  }
  std::vector<int> out;for(int j=0;j<k;++j)if(state[j]==2)out.push_back(j+1);return wrap(out);
}

// [[Rcpp::export]]
List spolo_extract_cpp(NumericMatrix spectra,NumericVector area,NumericVector center,
                       double fraction=.1,double pure_radius=.005,int nthreads=1) {
  const int m=spectra.nrow(),k=spectra.ncol();
  if(area.size()!=k||center.size()!=m||nthreads<1||fraction<=0||fraction>1)stop("Invalid spectral extraction controls.");
  NumericMatrix E(m,k);NumericVector mass(k),distance(k);LogicalVector pure(k);
  std::vector<int> ions;double weight=0;
  for(int c=0;c<m;++c)if(center[c]>0){ions.push_back(c);weight+=center[c];}
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(nthreads)
  #endif
  for(int j=0;j<k;++j) {
    if(area[j]<=0){distance[j]=R_PosInf;continue;}
    const double* s=&spectra(0,j);double* e=&E(0,j);
    distance[j]=js_distance(s,area[j],center.begin(),m);
    if(distance[j]<=pure_radius){pure[j]=true;for(int c=0;c<m;++c)e[c]=s[c];mass[j]=area[j];continue;}
    std::vector<std::pair<double,double> > ratio;ratio.reserve(ions.size());
    for(int c:ions)ratio.emplace_back(s[c]/center[c],center[c]);
    std::sort(ratio.begin(),ratio.end(),[](const std::pair<double,double>& a,const std::pair<double,double>& b){return a.first<b.first;});
    double cumulative=0,amplitude=0;
    for(const auto& item:ratio){cumulative+=item.second;amplitude=item.first;if(cumulative>=fraction*weight)break;}
    double total=0;for(int c:ions){e[c]=std::min(s[c],amplitude*center[c]);total+=e[c];}mass[j]=total;
  }
  return List::create(_["E"]=E,_["mass"]=mass,_["pure"]=pure,_["JS"]=distance);
}
