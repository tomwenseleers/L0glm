#include <Rcpp.h>
#include <algorithm>
#include <cmath>
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
        const double wa=cluster_weight(W,weights.size(),n,i,j), wb=cluster_weight(W,weights.size(),n,i,k);
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
