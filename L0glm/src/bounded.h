#ifndef L0GLM_BOUNDED_H
#define L0GLM_BOUNDED_H
#include <type_traits>
#include <vector>
#include <unordered_set>
#include <limits>
#include <cmath>
#include <algorithm>
namespace L0glm {
using Eigen::MatrixXd;using Eigen::SparseMatrix;using Eigen::VectorXd;
struct BoundedResult { VectorXd coefficients;bool success;int iterations;double error; };
template<class Matrix>
BoundedResult bounded_polish(Matrix G, VectorXd h, VectorXd x,
                         VectorXd lower, VectorXd upper, void (*interrupt)()=nullptr) {
  const int p = h.size();
  if (G.rows()!=p || G.cols()!=p || x.size()!=p || lower.size()!=p || upper.size()!=p)
    return BoundedResult{x,false,0,std::numeric_limits<double>::infinity()};
  const VectorXd original = x;
  VectorXd scale(p);
  std::vector<int> state(p), free;
  bool valid = h.allFinite() && x.allFinite();
  for (int j=0; j<p; ++j) {
    const double diagonal = G.coeff(j,j);
    valid = valid && std::isfinite(diagonal) && diagonal>=0 &&
      !std::isnan(lower[j]) && !std::isnan(upper[j]) && lower[j]<=upper[j];
    scale[j] = diagonal>0 ? std::sqrt(diagonal) : 1;
    x[j] = std::max(lower[j],std::min(upper[j],x[j]));
    // A zero-information Gram column has no identifiable update.
    if (diagonal==0) {
      valid = valid && h[j]==0;
      lower[j] = upper[j] = x[j];
    }
  }
  int iterations = 0;
  double error = std::numeric_limits<double>::infinity();
  bool success = valid && p==0;
  if (valid && p) {
    if constexpr (std::is_same_v<Matrix,MatrixXd>) {
      for (int j=0; j<p; ++j) for (int i=0; i<p; ++i)
        G(i,j) /= scale[i]*scale[j];
    } else {
      for (int j=0; j<G.outerSize(); ++j)
        for (typename Matrix::InnerIterator entry(G,j); entry; ++entry)
          entry.valueRef() /= scale[entry.row()]*scale[entry.col()];
    }
    h.array() /= scale.array();
    x.array() *= scale.array();
    lower.array() *= scale.array();
    upper.array() *= scale.array();
    const double denominator = std::max(1.,h.cwiseAbs().maxCoeff());
    const double tolerance = 1e-10*denominator;
    for (int j=0; j<p; ++j)
      state[j] = x[j]<=lower[j] ? -1 : x[j]>=upper[j] ? 1 : 0;

    for (; iterations<5*p+20 && x.allFinite(); ++iterations) {
      if (interrupt && iterations%64==0) interrupt();
      free.clear();
      for (int j=0; j<p; ++j) if (state[j]==0) free.push_back(j);
      VectorXd fixed = x;
      for (int j:free) fixed[j] = 0;
      const VectorXd rhs = h-G*fixed;
      const int q = free.size();
      VectorXd target(q);
      for (int k=0; k<q; ++k) target[k] = rhs[free[k]];
      if (q) {
        if constexpr (std::is_same_v<Matrix,MatrixXd>) {
          MatrixXd A(q,q);
          for (int j=0; j<q; ++j) for (int i=0; i<q; ++i)
            A(i,j) = G(free[i],free[j]);
          Eigen::LLT<MatrixXd> factor(A);
          if (factor.info()!=Eigen::Success) break;
          target = factor.solve(target).eval();
        } else {
          std::vector<int> index(p,-1);
          for (int k=0; k<q; ++k) index[free[k]] = k;
          std::vector<Eigen::Triplet<double>> entries;
          for (int k=0; k<q; ++k)
            for (typename Matrix::InnerIterator entry(G,free[k]); entry; ++entry)
              if (index[entry.row()]>=0)
                entries.emplace_back(index[entry.row()],k,entry.value());
          SparseMatrix<double> A(q,q);
          A.setFromTriplets(entries.begin(),entries.end());
          Eigen::SimplicialLLT<SparseMatrix<double>> factor(A);
          if (factor.info()!=Eigen::Success) break;
          target = factor.solve(target).eval();
        }
        if (!target.allFinite()) break;
        // Stop at the first bound, retain feasibility, then solve the new face.
        double step = 1;
        int hit = -1, side = 0;
        for (int k=0; k<q; ++k) {
          const int j = free[k];
          double fraction = 1;
          int bound = 0;
          if (target[k]<lower[j]) {
            fraction = (lower[j]-x[j])/(target[k]-x[j]); bound = -1;
          } else if (target[k]>upper[j]) {
            fraction = (upper[j]-x[j])/(target[k]-x[j]); bound = 1;
          }
          if (bound && fraction<=step) {
            step = std::max(0.,fraction); hit = j; side = bound;
          }
        }
        for (int k=0; k<q; ++k) x[free[k]] += step*(target[k]-x[free[k]]);
        if (hit>=0) {
          state[hit] = side;
          x[hit] = side<0 ? lower[hit] : upper[hit];
          continue;
        }
      }
      const VectorXd gradient = G*x-h;
      double violation = tolerance;
      int release = -1;
      for (int j=0; j<p; ++j) if (lower[j]<upper[j]) {
        const double value = state[j]<0 ? -gradient[j] : state[j]>0 ? gradient[j] : 0;
        if (value>violation) { violation = value; release = j; }
      }
      if (release<0) { success = true; break; }
      state[release] = 0;
    }
    // Fresh KKT residual, including stationarity on the free coefficients.
    for (int j=0; j<p; ++j) x[j] = std::max(lower[j],std::min(upper[j],x[j]));
    const VectorXd gradient = G*x-h;
    error = 0;
    for (int j=0; j<p; ++j) if (lower[j]<upper[j]) {
      double value = gradient[j];
      if ((x[j]<=lower[j] && value>0) || (x[j]>=upper[j] && value<0)) value = 0;
      error = std::max(error,std::abs(value));
    }
    error /= denominator;
    success = success && x.allFinite() && std::isfinite(error) && error<=1e-10;
    x.array() /= scale.array();
  }
  return BoundedResult{success ? x : original,success,iterations,p ? error : 0};
}

template<class Matrix> BoundedResult bounded_block(const Matrix& original_G,const VectorXd& original_h,
  const VectorXd& original_x,const VectorXd& original_lower,const VectorXd& original_upper,
  void (*interrupt)()=nullptr) {
  Matrix G=original_G;VectorXd h=original_h,x=original_x,lower=original_lower,upper=original_upper;
  const int p=h.size();
  if(G.rows()!=p||G.cols()!=p||x.size()!=p||lower.size()!=p||upper.size()!=p)
    return BoundedResult{original_x,false,0,std::numeric_limits<double>::infinity()};
  if(!p)return BoundedResult{x,true,0,0};
  VectorXd scale(p);std::vector<int> state(p),free;
  bool valid=h.allFinite()&&x.allFinite();
  for(int j=0;j<p;++j) {
    double d=G.coeff(j,j);
    valid=valid&&std::isfinite(d)&&d>=0&&!std::isnan(lower[j])&&!std::isnan(upper[j])&&lower[j]<=upper[j];
    scale[j]=d>0?std::sqrt(d):1;x[j]=std::max(lower[j],std::min(upper[j],x[j]));
    if(d==0){valid=valid&&h[j]==0;lower[j]=upper[j]=x[j];}
  }
  if(!valid)return bounded_polish(original_G,original_h,original_x,original_lower,original_upper,interrupt);
  if constexpr(std::is_same_v<Matrix,MatrixXd>) {
    for(int j=0;j<p;++j)for(int i=0;i<p;++i)G(i,j)/=scale[i]*scale[j];
  }else for(int j=0;j<G.outerSize();++j)for(typename Matrix::InnerIterator it(G,j);it;++it)
    it.valueRef()/=scale[it.row()]*scale[it.col()];
  h.array()/=scale.array();x.array()*=scale.array();lower.array()*=scale.array();upper.array()*=scale.array();
  const double denominator=std::max(1.,h.cwiseAbs().maxCoeff()),tol=1e-10*denominator;
  for(int j=0;j<p;++j)state[j]=x[j]<=lower[j]?-1:x[j]>=upper[j]?1:0;
  std::unordered_set<size_t> visited;
  for(int iteration=0;iteration<std::min(40,std::max(20,2*p));++iteration) {
    if(interrupt)interrupt();
    size_t hash=1469598103934665603ULL;
    for(int j=0;j<p;++j){hash^=(size_t)(state[j]+2);hash*=1099511628211ULL;}
    if(!visited.insert(hash).second)break;
    free.clear();for(int j=0;j<p;++j)if(state[j]==0)free.push_back(j);
    VectorXd fixed=x;for(int j:free)fixed[j]=0;
    VectorXd rhs=h-G*fixed;const int q=free.size();VectorXd target(q);
    for(int j=0;j<q;++j)target[j]=rhs[free[j]];
    if(q) {
      if constexpr(std::is_same_v<Matrix,MatrixXd>) {
        MatrixXd A(q,q);for(int j=0;j<q;++j)for(int i=0;i<q;++i)A(i,j)=G(free[i],free[j]);
        Eigen::LLT<MatrixXd> factor(A);if(factor.info()!=Eigen::Success)break;
        target=factor.solve(target).eval();
      }else {
        std::vector<int> index(p,-1);for(int k=0;k<q;++k)index[free[k]]=k;
        std::vector<Eigen::Triplet<double>> entries;
        for(int j=0;j<q;++j)for(typename Matrix::InnerIterator it(G,free[j]);it;++it)
          if(index[it.row()]>=0)entries.emplace_back(index[it.row()],j,it.value());
        SparseMatrix<double> A(q,q);A.setFromTriplets(entries.begin(),entries.end());
        Eigen::SimplicialLLT<SparseMatrix<double>> factor(A);if(factor.info()!=Eigen::Success)break;
        target=factor.solve(target).eval();
      }
      if(!target.allFinite())break;
    }
    bool changed=false;
    for(int k=0;k<q;++k) {
      const int j=free[k];x[j]=target[k];
      if(x[j]<lower[j]){state[j]=-1;x[j]=lower[j];changed=true;}
      else if(x[j]>upper[j]){state[j]=1;x[j]=upper[j];changed=true;}
    }
    if(changed)continue;
    VectorXd gradient=G*x-h;
    for(int j=0;j<p;++j)if(lower[j]<upper[j]) {
      if((state[j]<0&&gradient[j]<-tol)||(state[j]>0&&gradient[j]>tol)){state[j]=0;changed=true;}
    }
    if(changed)continue;
    double error=0;
    for(int j=0;j<p;++j)if(lower[j]<upper[j]) {
      double g=gradient[j];
      if((x[j]<=lower[j]&&g>0)||(x[j]>=upper[j]&&g<0))g=0;
      error=std::max(error,std::abs(g));
    }
    error/=denominator;
    if(x.allFinite()&&std::isfinite(error)&&error<=1e-10) {
      x.array()/=scale.array();return BoundedResult{x,true,iteration+1,error};
    }
    break;
  }
  VectorXd warm=x.cwiseMax(lower).cwiseMin(upper);
  warm.array()/=scale.array();
  return bounded_polish(original_G,original_h,warm,original_lower,original_upper,interrupt);
}

// Compare quadratics without subtracting a potentially large constant y'y.
template<class Matrix> bool bounded_improves(const Matrix& G,const VectorXd& h,
  const VectorXd& before,const VectorXd& after) {
  const VectorXd delta=after-before,gradient=G*before-h,curvature=G*delta;
  const double change=2*delta.dot(gradient)+delta.dot(curvature);
  const double magnitude=2*delta.cwiseAbs().dot(gradient.cwiseAbs())+
    delta.cwiseAbs().dot(curvature.cwiseAbs());
  return std::isfinite(change)&&change<=128*std::numeric_limits<double>::epsilon()*std::max(1.,magnitude);
}
inline bool nonnegative_polish_needed(const VectorXd& b,const VectorXd& lower,const VectorXd& upper) {
  return b.size()>0&&b.allFinite()&&(lower.array()>=0).any()&&
    ((b.array()<lower.array()).any()||(b.array()>upper.array()).any());
}
// Only callers that already own a predictor Gram matrix enter this function.
// BLAS LLT paths can supply its lower triangle; fill the other half on demand.
template<class Matrix> VectorXd nonnegative_polish(Matrix& G,const VectorXd& h,
  const VectorXd& b,const VectorXd& lower,const VectorXd& upper,bool lower_triangle=false) {
  if(!nonnegative_polish_needed(b,lower,upper))return b;
  if constexpr(std::is_same_v<Matrix,MatrixXd>) {
    if(lower_triangle)for(int j=0;j<G.cols();++j)for(int i=j+1;i<G.rows();++i)G(j,i)=G(i,j);
  }
  const VectorXd clipped=b.cwiseMax(lower).cwiseMin(upper);
  const auto answer=bounded_block(G,h,clipped,lower,upper);
  return answer.success&&bounded_improves(G,h,clipped,answer.coefficients)?answer.coefficients:clipped;
}
template<class Matrix> void nonnegative_polish_selected(Matrix& G,const VectorXd& h,VectorXd& b,
  const std::vector<int>& active,const VectorXd& scale,const VectorXd& lo,const VectorXd& hi,
  bool lower_triangle=false) {
  VectorXd lower(b.size()),upper(b.size());
  for(int j=0;j<b.size();++j) {
    const int a=active[j];
    lower[j]=lo[a]<=-std::numeric_limits<double>::max()?-std::numeric_limits<double>::infinity():lo[a]*scale[a];
    upper[j]=hi[a]>=std::numeric_limits<double>::max()?std::numeric_limits<double>::infinity():hi[a]*scale[a];
  }
  b=nonnegative_polish(G,h,b,lower,upper,lower_triangle);
}
}
#endif
