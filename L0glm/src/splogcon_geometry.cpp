#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

namespace {
double js_distance(const double* spectrum,double area,const double* p,int m) {
  double value=0;
  for(int c=0;c<m;++c) {
    const double q=spectrum[c]/area,mid=.5*(q+p[c]);
    if(q>0)value+=.5*q*std::log(q/mid);if(p[c]>0)value+=.5*p[c]*std::log(p[c]/mid);
  }
  return value;
}
}

// Linear-time least concave majorant shared by local and whole-design updates.
namespace {
double project_logconcave(const double* profile,double* out,int n,double min_relative) {
  double maximum=0;for(int i=0;i<n;++i)maximum=std::max(maximum,profile[i]);
  long double total=0;
  for(int i=0;i<n;++i){out[i]=std::max(profile[i],maximum*min_relative);total+=out[i];}
  for(int i=0;i<n;++i)out[i]/=(double)total;
  if(n<3)return 0;
  std::vector<int> stack;std::vector<double> y(n);
  for(int i=0;i<n;++i)y[i]=std::log(out[i]);
  for(int i=0;i<n;++i) {
    stack.push_back(i);
    while(stack.size()>=3) {
      const int k=stack.size(),a=stack[k-3],b=stack[k-2],c=stack[k-1];
      if((y[b]-y[a])/(b-a)>=(y[c]-y[b])/(c-b))break;
      stack[k-2]=c;stack.pop_back();
    }
  }
  total=0;
  for(size_t k=1;k<stack.size();++k) {
    const int a=stack[k-1],b=stack[k];
    for(int i=a;i<b;++i){out[i]=std::exp(y[a]+(y[b]-y[a])*(double)(i-a)/(b-a));total+=out[i];}
  }
  out[n-1]=std::exp(y[n-1]);total+=out[n-1];
  const double gap=std::max(0.,(double)total-1);
  for(int i=0;i<n;++i)out[i]/=(double)total;
  return gap;
}
}

// [[Rcpp::export]]
List slcn_logconcave_cpp(NumericVector profile,double min_relative=1e-12) {
  const int n=profile.size();
  if(!n||!std::isfinite(min_relative)||min_relative<0||min_relative>1)stop("Invalid log-concavity projection.");
  double maximum=0;for(double v:profile){if(!std::isfinite(v))stop("Profiles must be finite.");maximum=std::max(maximum,v);}
  if(maximum<=0)stop("A profile needs positive mass.");
  NumericVector out(n);
  const double gap=project_logconcave(profile.begin(),out.begin(),n,min_relative);
  return List::create(_["h"]=out,_["gap"]=gap);
}

// One copy of the sparse design; all column projections write disjoint slots.
// [[Rcpp::export]]
S4 slcn_project_design_cpp(S4 design,int nthreads=1,double min_relative=1e-12) {
  if(!Rf_inherits(design,"dgCMatrix")||nthreads<1||!std::isfinite(min_relative)||min_relative<0||min_relative>1)
    stop("Invalid sparse projection controls.");
  S4 out=clone(design);IntegerVector ptr=out.slot("p"),index=out.slot("i"),dims=out.slot("Dim");
  NumericVector values=out.slot("x");const int p=dims[1];
  for(double v:values)if(!std::isfinite(v)||v<0)stop("Profiles must be finite and nonnegative.");
  for(int j=0;j<p;++j) {
    if(ptr[j]==ptr[j+1])stop("A profile needs positive mass.");
    double maximum=0;
    for(int z=ptr[j];z<ptr[j+1];++z) {
      maximum=std::max(maximum,values[z]);
      if(z>ptr[j]&&index[z]!=index[z-1]+1)stop("Profiles need contiguous support.");
    }
    if(maximum<=0)stop("A profile needs positive mass.");
  }
  double* x=values.begin();
  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic,1) num_threads(nthreads)
  #endif
  for(int j=0;j<p;++j)project_logconcave(x+ptr[j],x+ptr[j],ptr[j+1]-ptr[j],min_relative);
  return out;
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
