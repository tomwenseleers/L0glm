#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::plugins(cpp17)]]

// Gibbs conditionals retain the complete precision matrix. Tail probabilities
// are sampled on the log-survival scale to avoid inverse-CDF underflow.
// [[Rcpp::export]]
NumericMatrix posterior_truncated_cpp(NumericVector mean, NumericMatrix precision,
    LogicalVector constrained, NumericVector start, int nsim, int burnin, int thin) {
  int k=mean.size();
  if(precision.nrow()!=k || precision.ncol()!=k || constrained.size()!=k ||
     start.size()!=k || nsim<1 || burnin<0 || thin<1)stop("Invalid posterior sampler dimensions.");
  NumericVector x=clone(start), residual(k);
  NumericMatrix out(k,nsim);
  for(int j=0;j<k;++j)if(!R_finite(precision(j,j)) || precision(j,j)<=0)
    stop("Posterior precision must have a positive finite diagonal.");
  int saved=0,total=burnin+nsim*thin;
  for(int iteration=0;iteration<total;++iteration) {
    if(iteration%32==0) {
      checkUserInterrupt();
      for(int i=0;i<k;++i) {
        residual[i]=0;
        for(int j=0;j<k;++j)residual[i]+=precision(i,j)*(x[j]-mean[j]);
      }
    }
    for(int j=0;j<k;++j) {
      double sd=1/std::sqrt(precision(j,j));
      double center=x[j]-residual[j]/precision(j,j);
      double value;
      if(constrained[j]) {
        double a=-center/sd;
        double log_tail=R::pnorm(a,0,1,false,true);
        double log_u=std::log(R::runif(0,1));
        value=center+sd*R::qnorm(log_tail+log_u,0,1,false,true);
        if(value<0 && value>-1e-12*sd)value=0; // roundoff at an extreme boundary
      } else value=R::rnorm(center,sd);
      if(!R_finite(value) || (constrained[j] && value<0))stop("Nonfinite truncated-normal draw.");
      double change=value-x[j];x[j]=value;
      for(int i=0;i<k;++i)residual[i]+=precision(i,j)*change;
    }
    if(iteration>=burnin && (iteration-burnin+1)%thin==0) {
      for(int j=0;j<k;++j)out(j,saved)=x[j];
      ++saved;
    }
  }
  return out;
}
