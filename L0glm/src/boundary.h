#ifndef L0GLM_BOUNDARY_H
#define L0GLM_BOUNDARY_H

namespace L0glm {
// Conditional endpoint tests share the same arithmetic for both storage types.
template<class Function>
inline void boundary_column(const Eigen::MappedSparseMatrix<double>& X,int j,Function f) {
  for(Eigen::MappedSparseMatrix<double>::InnerIterator it(X,j);it;++it)
    if(it.value()!=0.0) f(it.row(),it.value());
}
template<class Function>
inline void boundary_column(const Eigen::Map<Eigen::MatrixXd>& X,int j,Function f) {
  for(int i=0;i<X.rows();++i) if(X(i,j)!=0.0) f(i,X(i,j));
}

inline bool boundary_nonnegative(const Eigen::MappedSparseMatrix<double>& X) {
  for(int j=0;j<X.outerSize();++j)
    for(Eigen::MappedSparseMatrix<double>::InnerIterator it(X,j);it;++it)
      if(!(it.value()>=0.0)) return false;
  return true;
}
inline bool boundary_nonnegative(const Eigen::Map<Eigen::MatrixXd>& X) {
  return (X.array()>=0.0).all();
}

// family=0: Gaussian identity; family=1: Poisson identity.
// Only zero-bound positive coefficients are eligible. The candidate must be
// feasible and improve the unpenalized objective; no size threshold is used.
template<class Design>
int boundary_identity(const Design& X,const Eigen::VectorXd& y,
  const Eigen::VectorXd& prior,Eigen::VectorXd& beta,Eigen::VectorXd& eta,
  const Eigen::VectorXd& lower,int family,bool nonnegative_components=false) {
  int removed=0;
  for(int j=0;j<beta.size();++j) if(lower[j]==0.0 && beta[j]>0.0) {
    const double coefficient=beta[j];
    double score=0.0,change=0.0;bool feasible=true;
    boundary_column(X,j,[&](int i,double x) {
      if(!feasible) return;
      double next=eta[i]-x*coefficient;
      if(!std::isfinite(next)) {feasible=false;return;}
      if(family==1 && next<0.0) {
        // With nonnegative components the remainder is nonnegative exactly;
        // repeated subtraction can otherwise leave accumulated roundoff.
        if(!nonnegative_components) {feasible=false;return;}
        next=0.0;
      }
      const double w=prior[i];if(w==0.0) return;
      if(family==0) {
        score += w*x*(next-y[i]);
        change += w*(eta[i]-next)*(2*y[i]-eta[i]-next)/2.0;
      } else {
        if((y[i]>0.0 && next<=0.0) || eta[i]<0.0) {feasible=false;return;}
        score += w*x*(y[i]==0.0 ? 1.0 : 1.0-y[i]/next);
        change -= w*x*coefficient;
        if(y[i]>0.0) change -= w*y[i]*std::log1p(-x*coefficient/eta[i]);
      }
    });
    if(feasible && score>=0.0 && std::isfinite(change) && change<=0.0) {
      beta[j]=0.0;++removed;
      boundary_column(X,j,[&](int i,double x) {
        eta[i]-=x*coefficient;
        if(family==1) eta[i]=std::max(0.0,eta[i]);
      });
    }
  }
  return removed;
}
} // namespace L0glm
#endif
