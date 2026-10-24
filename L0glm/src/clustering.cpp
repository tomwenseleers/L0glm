#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <queue>
#include <vector>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

namespace {
inline double cluster_weight(const double* w, R_xlen_t length, int n, int i, int j) {
  return w[length == 1 ? 0 : length == n ? i : i + static_cast<R_xlen_t>(n)*j];
}
void cluster_dimensions(const NumericMatrix& y, const NumericVector& w, int threads) {
  const R_xlen_t size = static_cast<R_xlen_t>(y.nrow())*y.ncol();
  if (y.nrow()<1 || y.ncol()<1 || threads<1 ||
      !(w.size()==1 || w.size()==y.nrow() || w.size()==size))
    stop("Invalid clustering dimensions, weights or thread count.");
}
}

// Smoothed Poisson log-likelihood gain with fixed channel size factors.
// Workers write separate entries in R's packed lower-triangular dist layout.
// [[Rcpp::export]]
NumericVector L0glm_poissondist_cpp(NumericMatrix y, NumericVector weights,
                                   NumericVector size_fac, int nthreads=1) {
  cluster_dimensions(y,weights,nthreads);
  const int n=y.nrow(), m=y.ncol();
  const R_xlen_t nw=weights.size();
  if (size_fac.size()!=m) stop("One size factor is required per channel.");
  for (double v : y) if (!std::isfinite(v) || v<0) stop("Poisson observations must be finite and nonnegative.");
  for (double v : weights) if (!std::isfinite(v) || v<0) stop("Weights must be finite and nonnegative.");
  for (int j=0;j<m;++j) {
    if (!std::isfinite(size_fac[j]) || size_fac[j]<0) stop("Invalid size factors.");
    if (size_fac[j]==0) for (int i=0;i<n;++i)
      if (y(i,j)>0 && cluster_weight(weights.begin(),weights.size(),n,i,j)>0)
        stop("A positive-count channel needs a positive size factor.");
  }
  const R_xlen_t count=static_cast<R_xlen_t>(m)*(m-1)/2;
  NumericVector result(count);
  const double *Y=y.begin(), *W=weights.begin(), *S=size_fac.begin();
  double *out=result.begin();
  nthreads=std::min(nthreads,std::max(1,m-1));
#ifdef _OPENMP
  nthreads=std::min(nthreads,omp_get_num_procs());
#endif
#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
#endif
  for (int j=0;j<m-1;++j) {
    const R_xlen_t start=static_cast<R_xlen_t>(j)*(2.0*m-j-1)/2;
    for (int k=j+1;k<m;++k) {
      double total=0;
      for (int i=0;i<n;++i) {
        const double a=Y[i+static_cast<R_xlen_t>(n)*j], b=Y[i+static_cast<R_xlen_t>(n)*k];
        const double wa=cluster_weight(W,nw,n,i,j), wb=cluster_weight(W,nw,n,i,k);
        const double den=wa*S[j]+wb*S[k];
        if (!(den>0)) continue;
        const double centre=(wa*a+wb*b)/den, ma=S[j]*centre, mb=S[k]*centre;
        if (wa>0 && ma>0) {
          const double d=(wa*a+1)/(wa*ma+1);
          total+=wa*(ma-ma*d+(a>0 ? a*std::log(d) : 0));
        }
        if (wb>0 && mb>0) {
          const double d=(wb*b+1)/(wb*mb+1);
          total+=wb*(mb-mb*d+(b>0 ? b*std::log(d) : 0));
        }
      }
      out[start+k-j-1]=std::max(0.0,total);
    }
  }
  return result;
}

namespace {
// Unit deviances for independent scalar mean models. Dispersion is absorbed
// into precision weights before dispatch; NB shape stays fixed.
inline double cluster_deviance(double y, double mu, int family, double theta) {
  if (y==mu) return 0;
  if (family==2) return (y-mu)*(y-mu);
  double value=0;
  if (family==1) value=2*((y>0 ? y*std::log(y/mu) : 0)-y+mu);
  else if (family==3) value=2*((y>0 ? y*std::log(y/mu) : 0)+
    (y<1 ? (1-y)*std::log((1-y)/(1-mu)) : 0));
  else if (family==4) {const double r=(y-mu)/mu;value=2*(r-std::log1p(r));}
  else if (family==5) value=(y-mu)*(y-mu)/(y*mu*mu);
  else value=2*((y>0 ? y*std::log(y/mu) : 0)-(y+theta)*std::log((y+theta)/(mu+theta)));
  return std::max(0.0,value);
}
inline int cluster_threads(int requested, int jobs) {
  int value=std::min(requested,std::max(1,jobs));
#ifdef _OPENMP
  value=std::min(value,omp_get_num_procs());
#endif
  return value;
}
void cluster_domain(const NumericMatrix& y,int family,double theta) {
  if(family<1 || family>6 || (family==6 && (!std::isfinite(theta)||theta<=0)))
    stop("Unsupported scalar family or NB shape.");
  for(double v:y) if(!std::isfinite(v) || ((family==1||family==6)&&v<0) ||
    (family==3&&(v<0||v>1)) || ((family==4||family==5)&&v<=0))
    stop("Observations are outside the scalar family's domain.");
}
struct ClusterPair {
  double loss;int a,b;
  bool operator<(const ClusterPair& z) const {
    if(loss!=z.loss) return loss>z.loss;
    return a!=z.a ? a>z.a : b>z.b;
  }
};
}

// [[Rcpp::export]]
NumericVector L0glm_deviance_dist_cpp(NumericMatrix y,NumericVector weights,
                                      int family,double theta,int nthreads=1) {
  cluster_dimensions(y,weights,nthreads);cluster_domain(y,family,theta);
  const int n=y.nrow(),m=y.ncol();const R_xlen_t nw=weights.size();
  NumericVector result(static_cast<R_xlen_t>(m)*(m-1)/2);
  const double *Y=y.begin(),*W=weights.begin();double *out=result.begin();
  nthreads=cluster_threads(nthreads,m-1);
  for(int block=0;block<m-1;block+=32) {
    checkUserInterrupt();const int end=std::min(m-1,block+32);
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
#endif
    for(int j=block;j<end;++j) {
      const R_xlen_t start=static_cast<R_xlen_t>(j)*(2.0*m-j-1)/2;
      for(int k=j+1;k<m;++k) {
        double loss=0;
        for(int i=0;i<n;++i) {
          const double a=Y[i+static_cast<R_xlen_t>(n)*j],b=Y[i+static_cast<R_xlen_t>(n)*k];
          const double wa=cluster_weight(W,nw,n,i,j),wb=cluster_weight(W,nw,n,i,k),total=wa+wb;
          if(total<=0) continue;
          const double mu=(wa*a+wb*b)/total;
          if(wa>0) loss+=wa*cluster_deviance(a,mu,family,theta);
          if(wb>0) loss+=wb*cluster_deviance(b,mu,family,theta);
        }
        out[start+k-j-1]=loss;
      }
    }
  }
  return result;
}

// Exact agglomerative merge cost for the specified pooled-mean likelihood.
// For Poisson size factors, S=sum(w*y), T=sum(w*size), and S/T is a rate.
// [[Rcpp::export]]
List L0glm_bregman_tree_cpp(NumericMatrix y,NumericVector weights,
                            NumericVector sizes,int family,double theta,int nthreads=1) {
  cluster_dimensions(y,weights,nthreads);cluster_domain(y,family,theta);
  const int n=y.nrow(),m=y.ncol(),nodes=2*m-1;const R_xlen_t nw=weights.size();
  if(sizes.size()!=m) stop("Incompatible size factors.");
  const double *Y=y.begin(),*W=weights.begin(),*scale=sizes.begin();
  std::vector<double> sums(static_cast<R_xlen_t>(n)*nodes),mass(sums.size());
  std::vector<unsigned char> active(nodes,0);std::vector<int> first(nodes);
  for(int j=0;j<m;++j) {
    active[j]=1;first[j]=j;
    for(int i=0;i<n;++i) {
      const double w=cluster_weight(W,nw,n,i,j);const R_xlen_t pos=i+static_cast<R_xlen_t>(n)*j;
      sums[pos]=w*Y[pos];mass[pos]=w*scale[j];
    }
  }
  auto cost=[&](int a,int b) {
    double loss=0;
    for(int i=0;i<n;++i) {
      const R_xlen_t ia=i+static_cast<R_xlen_t>(n)*a,ib=i+static_cast<R_xlen_t>(n)*b;
      const double wa=mass[ia],wb=mass[ib],total=wa+wb;
      if(total<=0) continue;
      const double mu=(sums[ia]+sums[ib])/total;
      if(wa>0) loss+=wa*cluster_deviance(sums[ia]/wa,mu,family,theta);
      if(wb>0) loss+=wb*cluster_deviance(sums[ib]/wb,mu,family,theta);
    }
    return loss;
  };
  std::priority_queue<ClusterPair> queue;
  std::vector<double> losses(m);nthreads=cluster_threads(nthreads,m-1);
  for(int a=0;a<m-1;++a) {
    checkUserInterrupt();
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
    for(int b=a+1;b<m;++b) losses[b]=cost(a,b);
    for(int b=a+1;b<m;++b) queue.push({losses[b],a,b});
  }
  IntegerMatrix merge(std::max(0,m-1),2);NumericVector height(std::max(0,m-1));
  auto code=[&](int node){return node<m ? -node-1 : node-m+1;};
  for(int step=0;step<m-1;++step) {
    checkUserInterrupt();
    while(!queue.empty() && !(active[queue.top().a]&&active[queue.top().b])) queue.pop();
    if(queue.empty() || !std::isfinite(queue.top().loss)) stop("Nonfinite clustering loss; rescale the observations.");
    const ClusterPair pair=queue.top();queue.pop();int a=pair.a,b=pair.b,newnode=m+step;
    if(first[a]>first[b]) std::swap(a,b);
    merge(step,0)=code(a);merge(step,1)=code(b);height[step]=pair.loss;
    active[a]=active[b]=0;active[newnode]=1;first[newnode]=std::min(first[a],first[b]);
    for(int i=0;i<n;++i) {
      const R_xlen_t dest=i+static_cast<R_xlen_t>(n)*newnode;
      sums[dest]=sums[i+static_cast<R_xlen_t>(n)*a]+sums[i+static_cast<R_xlen_t>(n)*b];
      mass[dest]=mass[i+static_cast<R_xlen_t>(n)*a]+mass[i+static_cast<R_xlen_t>(n)*b];
    }
    std::vector<int> others;for(int j=0;j<newnode;++j) if(active[j]) others.push_back(j);
    std::vector<double> next(others.size());
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
    for(int t=0;t<static_cast<int>(others.size());++t) next[t]=cost(others[t],newnode);
    for(int t=0;t<static_cast<int>(others.size());++t) queue.push({next[t],others[t],newnode});
  }
  std::vector<int> order,stack(1,nodes-1);
  while(!stack.empty()) {
    const int node=stack.back();stack.pop_back();
    if(node<m) order.push_back(node+1);
    else {const int row=node-m;for(int side=1;side>=0;--side) {
      const int value=merge(row,side);stack.push_back(value<0 ? -value-1 : m+value-1);
    }}
  }
  return List::create(_["merge"]=merge,_["height"]=height,_["order"]=wrap(order));
}

// Profile deviance for candidate cuts, always evaluated on the same channels.
// [[Rcpp::export]]
List L0glm_cluster_profile_cpp(NumericMatrix y,NumericVector weights,
                                NumericVector sizes,IntegerMatrix labels,
                                int family,double theta,int nthreads=1) {
  cluster_dimensions(y,weights,nthreads);cluster_domain(y,family,theta);
  const int n=y.nrow(),m=y.ncol(),q=labels.ncol();const R_xlen_t nw=weights.size();
  if(labels.nrow()!=m || sizes.size()!=m || q<1) stop("Invalid tree-cut labels.");
  std::vector<int> counts(q,0);
  for(int t=0;t<q;++t) for(int j=0;j<m;++j) {
    if(labels(j,t)<1 || labels(j,t)>m) stop("Invalid cluster labels.");
    counts[t]=std::max(counts[t],labels(j,t));
  }
  NumericVector losses(q),df(q);NumericMatrix centre(q==1 ? n : 0,q==1 ? counts[0] : 0);
  const double *Y=y.begin(),*W=weights.begin(),*scale=sizes.begin();const int *lab=labels.begin();
  double *out=losses.begin(),*degrees=df.begin(),*C=centre.begin();nthreads=cluster_threads(nthreads,q);
#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
#endif
  for(int t=0;t<q;++t) {
    const int k=counts[t];std::vector<double> sums(static_cast<R_xlen_t>(n)*k),mass(sums.size());
    for(int j=0;j<m;++j) for(int i=0;i<n;++i) {
      const int g=lab[j+m*t]-1;const R_xlen_t pos=i+static_cast<R_xlen_t>(n)*g;
      const double w=cluster_weight(W,nw,n,i,j);
      sums[pos]+=w*Y[i+static_cast<R_xlen_t>(n)*j];mass[pos]+=w*scale[j];
    }
    for(R_xlen_t z=0;z<static_cast<R_xlen_t>(n)*k;++z) if(mass[z]>0) {sums[z]/=mass[z];degrees[t]+=1;}
    for(int j=0;j<m;++j) for(int i=0;i<n;++i) {
      const double w=cluster_weight(W,nw,n,i,j);if(w<=0) continue;
      const double mu=scale[j]*sums[i+static_cast<R_xlen_t>(n)*(lab[j+m*t]-1)];
      out[t]+=w*cluster_deviance(Y[i+static_cast<R_xlen_t>(n)*j],mu,family,theta);
    }
    if(q==1) std::copy(sums.begin(),sums.end(),C);
  }
  return List::create(_["deviance"]=losses,_["df"]=df,_["centres"]=centre);
}
