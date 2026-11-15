#ifndef L0GLM_MEMORY_BUDGET_H
#define L0GLM_MEMORY_BUDGET_H
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace L0glm {
double available_memory();
inline double workspace_budget(double requested=-1,int threads=1) {
  return (requested<0?.25*available_memory():requested)/std::max(1,threads);
}
template<class Matrix> double matrix_workspace(const Matrix& A) {
  if constexpr(std::is_same_v<typename Matrix::StorageKind,Eigen::Sparse>)
    return 6*(12.*A.nonZeros()+4.*(A.cols()+1));
  else return 48.*A.rows()*A.cols();
}
template<class Matrix> double gram_workspace(const Matrix& A,bool rows=false) {
  const double dimension=rows?A.rows():A.cols();
  if constexpr(std::is_same_v<typename Matrix::StorageKind,Eigen::Sparse>) {
    double products=0;
    if(rows)for(int j=0;j<A.outerSize();++j) {
      double count=0;for(typename Matrix::InnerIterator it(A,j);it;++it)count+=it.value()!=0;
      products+=count*count;
    }else {
      std::vector<int> counts(A.rows(),0);
      for(int j=0;j<A.outerSize();++j)for(typename Matrix::InnerIterator it(A,j);it;++it)
        if(it.value()!=0)++counts[it.row()];
      for(double count:counts)products+=count*count;
    }
    return 72.*std::min(dimension*dimension,products)+24.*(dimension+1);
  }else return 48.*dimension*dimension;
}
inline void require_workspace(double estimate,double budget,const char* operation) {
  if(estimate>budget)throw std::runtime_error(std::string(operation)+
    " exceeds the Gram/factor workspace budget; use screening, Eigen CG/LSCG with diagonal preconditioning, or increase controlSolverGen(max.memory=...).");
}
struct MatrixFreeResult {
  Eigen::MatrixXd coefficients;
  bool success=true;int iterations=0;double error=0;
};
// Symmetric scaling and an independently recomputed true residual certify each RHS.
template<class Operator> MatrixFreeResult matrixfree_spd(Operator apply,
    const Eigen::MatrixXd& rhs,const Eigen::VectorXd& diagonal,const Eigen::MatrixXd& start,
    double tol=1e-10,int maxit=0) {
  const int n=rhs.rows();if(maxit<=0)maxit=std::max(100,2*n);
  if(!(tol>0))tol=1e-10;
  Eigen::VectorXd scale(n);for(int i=0;i<n;++i)scale[i]=diagonal[i]>0?std::sqrt(diagonal[i]):1;
  auto scaled=[&](const Eigen::VectorXd& v){return Eigen::VectorXd(apply(Eigen::VectorXd(v.array()/scale.array())).array()/scale.array());};
  MatrixFreeResult out;out.coefficients=Eigen::MatrixXd::Zero(n,rhs.cols());
  for(int j=0;j<rhs.cols();++j) {
    Eigen::VectorXd b=rhs.col(j).array()/scale.array(),x=start.col(j).array()*scale.array();
    Eigen::VectorXd r=b-scaled(x);if(r.squaredNorm()>b.squaredNorm()){x.setZero();r=b;}
    Eigen::VectorXd direction=r;double rr=r.squaredNorm(),error=0;
    const double norm_rhs=rhs.col(j).norm(),norm_b=b.norm();bool success=false;int iteration=0;
    if(norm_rhs==0){x.setZero();success=true;}
    else for(;iteration<maxit;++iteration) {
      if(std::sqrt(rr)<=tol*std::max(norm_b,std::numeric_limits<double>::min())) {
        Eigen::VectorXd value=x.array()/scale.array(),fresh=rhs.col(j)-apply(value);
        error=fresh.norm()/norm_rhs;
        if(std::isfinite(error)&&error<=tol){success=true;break;}
        r=fresh.array()/scale.array();direction=r;rr=r.squaredNorm();
      }
      Eigen::VectorXd product=scaled(direction);double curvature=direction.dot(product);
      if(!(curvature>0)||!std::isfinite(curvature))break;
      double alpha=rr/curvature;x+=alpha*direction;r-=alpha*product;
      double next=r.squaredNorm();direction=r+(next/rr)*direction;rr=next;
      if(!x.allFinite()||!std::isfinite(rr))break;
    }
    out.coefficients.col(j)=x.array()/scale.array();out.iterations+=iteration;
    if(norm_rhs>0)error=(rhs.col(j)-apply(Eigen::VectorXd(out.coefficients.col(j)))).norm()/norm_rhs;
    success=out.coefficients.col(j).allFinite()&&std::isfinite(error)&&error<=tol;
    out.success=out.success&&success;out.error=std::max(out.error,error);
  }
  return out;
}
template<class Matrix> MatrixFreeResult matrixfree_ridge(const Matrix& A,const Eigen::VectorXd& y,
    const Eigen::VectorXd& d,const Eigen::VectorXd& center,const Eigen::VectorXd& start,
    double tol=1e-10,int maxit=0) {
  const int p=A.cols();Eigen::VectorXd diagonal=d;
  if constexpr(std::is_same_v<typename Matrix::StorageKind,Eigen::Sparse>) {
    for(int j=0;j<A.outerSize();++j)for(typename Matrix::InnerIterator it(A,j);it;++it)
      diagonal[j]+=it.value()*it.value();
  }else for(int j=0;j<p;++j)diagonal[j]+=A.col(j).squaredNorm();
  Eigen::VectorXd rhs=A.transpose()*y;rhs+=d.cwiseProduct(center);
  auto apply=[&](const Eigen::VectorXd& v){return Eigen::VectorXd(A.transpose()*(A*v)+d.cwiseProduct(v));};
  return matrixfree_spd(apply,Eigen::MatrixXd(rhs),diagonal,Eigen::MatrixXd(start),tol,maxit);
}
// Observation-space solve, including the Schur complement of unpenalized terms.
template<class Matrix> MatrixFreeResult matrixfree_dual(const Matrix& A,const Eigen::VectorXd& y,
    const Eigen::VectorXd& d,Eigen::VectorXd center,double tol=1e-10,int maxit=0) {
  const int n=A.rows(),p=A.cols();const bool penalized=(d.array()>0).any();
  Eigen::VectorXd inverse=Eigen::VectorXd::Ones(p),diagonal=Eigen::VectorXd::Zero(n);
  std::vector<int> free;
  for(int j=0;j<p;++j)if(penalized) {
    inverse[j]=d[j]>0?1/d[j]:0;
    if(d[j]==0){free.push_back(j);center[j]=0;}
  }else center[j]=0;
  if constexpr(std::is_same_v<typename Matrix::StorageKind,Eigen::Sparse>) {
    for(int j=0;j<A.outerSize();++j)for(typename Matrix::InnerIterator it(A,j);it;++it)
      diagonal[it.row()]+=inverse[j]*it.value()*it.value();
  }else for(int j=0;j<p;++j)diagonal.array()+=inverse[j]*A.col(j).array().square();
  if(penalized)diagonal.array()+=1;
  auto apply=[&](const Eigen::VectorXd& v){
    Eigen::VectorXd value=A*(inverse.array()*(A.transpose()*v).array()).matrix();
    if(penalized)value+=v;
    return value;
  };
  Eigen::MatrixXd rhs(n,free.size()+1),U(n,free.size());rhs.col(0)=y-A*center;
  for(int j=0;j<(int)free.size();++j){U.col(j)=A.col(free[j]);rhs.col(j+1)=U.col(j);}
  auto out=matrixfree_spd(apply,rhs,diagonal,Eigen::MatrixXd::Zero(n,rhs.cols()),tol,maxit);
  Eigen::VectorXd u=out.coefficients.col(0),alpha(free.size());
  if(!free.empty()) {
    Eigen::MatrixXd T=out.coefficients.rightCols(free.size()),G=U.transpose()*T;
    alpha=G.colPivHouseholderQr().solve(U.transpose()*u);u.noalias()-=T*alpha;
  }
  out.coefficients=center+inverse.cwiseProduct(A.transpose()*u);
  for(int j=0;j<(int)free.size();++j)out.coefficients(free[j],0)=alpha[j];
  Eigen::VectorXd beta=out.coefficients.col(0),rhs_coefficient=A.transpose()*y+d.cwiseProduct(center);
  Eigen::VectorXd gradient=A.transpose()*(A*beta-y)+d.cwiseProduct(beta-center);
  const double residual=gradient.norm()/std::max(rhs_coefficient.norm(),std::numeric_limits<double>::min());
  if(!out.success||!std::isfinite(residual)||residual>tol) {
    auto refined=matrixfree_ridge(A,y,d,center,beta,tol,maxit);
    refined.iterations+=out.iterations;return refined;
  }
  out.error=residual;out.success=out.coefficients.allFinite();return out;
}
}
#endif
