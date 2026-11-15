#ifndef L0GLM_RIDGE_WORKSPACE_H
#define L0GLM_RIDGE_WORKSPACE_H
#include <memory>
#include <cstring>
#include <type_traits>
#include "memory_budget.h"

namespace L0glm {
// One workspace per active channel pass. No R objects or global mutable state.
template<class Matrix> struct RidgeWorkspace {
  Matrix A,base_gram;
  Eigen::VectorXd response,base_rhs;double gram_estimate=-1;
  std::vector<int> rows,columns;
  Eigen::VectorXd scales,last_weights,last_response,unweighted;
  bool ready=false,gram_ready=false,rhs_ready=false,lower_gram=false,inverse_scale=false,natural_mode=false;
  std::vector<int> factor_inner,factor_outer;
  using AMD=Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Lower,Eigen::AMDOrdering<int>>;
  using Natural=Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Lower,Eigen::NaturalOrdering<int>>;
  std::unique_ptr<AMD> amd;
  std::unique_ptr<Natural> natural;

  template<class Design> bool update(const Design& X,const std::vector<int>& selected_rows,
      const std::vector<int>& selected_columns,const Eigen::VectorXd& weights,
      const Eigen::VectorXd& scale,const Eigen::VectorXd& z,bool inverse=false) {
    const bool same=ready&&rows==selected_rows&&columns==selected_columns&&
      scales.size()==scale.size()&&(scales.array()==scale.array()).all()&&inverse_scale==inverse;
    const bool same_weights=same&&last_weights.size()==weights.size()&&
      (last_weights.array()==weights.array()).all();
    const bool same_response=same&&last_response.size()==z.size()&&
      (last_response.array()==z.array()).all();
    if(!same) {
      gram_estimate=-1;rows=selected_rows;columns=selected_columns;scales=scale;inverse_scale=inverse;
      if constexpr(std::is_same_v<Matrix,Eigen::SparseMatrix<double>>) {
        std::vector<int> map(X.rows(),-1);
        for(int i=0;i<(int)rows.size();++i)map[rows[i]]=i;
        std::vector<Eigen::Triplet<double>> entries;
        for(int j=0;j<(int)columns.size();++j)
          for(typename Design::InnerIterator it(X,columns[j]);it;++it)
            if(map[it.row()]>=0&&it.value()!=0)entries.emplace_back(map[it.row()],j,it.value());
        A.resize(rows.size(),columns.size());A.setFromTriplets(entries.begin(),entries.end());
        A.makeCompressed();unweighted=Eigen::Map<Eigen::VectorXd>(A.valuePtr(),A.nonZeros());
      }else A.resize(rows.size(),columns.size());
      factor_inner.clear();factor_outer.clear();
    }
    if(!same_weights) {
      Eigen::VectorXd root(weights.size());
      for(int i=0;i<weights.size();++i)root[i]=weights[i]>0?std::sqrt(weights[i]):0;
      if constexpr(std::is_same_v<Matrix,Eigen::SparseMatrix<double>>) {
        for(int j=0;j<A.cols();++j)for(int at=A.outerIndexPtr()[j];at<A.outerIndexPtr()[j+1];++at) {
          double value=unweighted[at]*root[rows[A.innerIndexPtr()[at]]];
          A.valuePtr()[at]=inverse?value*(1/scales[j]):value/scales[j];
        }
      }else for(int j=0;j<A.cols();++j)for(int i=0;i<A.rows();++i) {
        double value=X(rows[i],columns[j])*root[rows[i]];
        A(i,j)=inverse?value*(1/scales[j]):value/scales[j];
      }
      gram_ready=rhs_ready=false;last_weights=weights;
    }
    if(!same_weights||!same_response) {
      response.resize(rows.size());
      for(int i=0;i<(int)rows.size();++i)response[i]=z[rows[i]]*std::sqrt(weights[rows[i]]);
      rhs_ready=false;last_response=z;
    }
    ready=true;return !rows.empty();
  }
  double estimate() {if(gram_estimate<0)gram_estimate=gram_workspace(A);return gram_estimate;}
  const Matrix& gram(bool lower=false) {
    if(!gram_ready||lower_gram!=lower) {
      if constexpr(std::is_same_v<Matrix,Eigen::MatrixXd>) {
        const int n=A.rows(),p=A.cols();base_gram=Matrix::Zero(p,p);double one=1,zero=0;
        if(n&&p)F77_CALL(dsyrk)("L","T",&p,&n,&one,A.data(),&n,&zero,base_gram.data(),&p FCONE FCONE);
        if(!lower)base_gram.template triangularView<Eigen::Upper>()=
          base_gram.transpose().template triangularView<Eigen::Upper>();
      }else base_gram=A.transpose()*A;
      gram_ready=true;lower_gram=lower;
    }
    return base_gram;
  }
  const Eigen::VectorXd& rhs(bool gemv=false) {
    if(!rhs_ready) {
      if constexpr(std::is_same_v<Matrix,Eigen::MatrixXd>) {
        const int n=A.rows(),p=A.cols(),m=1,inc=1;double one=1,zero=0;base_rhs=Eigen::VectorXd::Zero(p);
        if(n&&p) {
          if(gemv)F77_CALL(dgemv)("T",&n,&p,&one,A.data(),&n,response.data(),&inc,&zero,base_rhs.data(),&inc FCONE);
          else F77_CALL(dgemm)("T","N",&p,&m,&n,&one,A.data(),&n,response.data(),&n,&zero,base_rhs.data(),&p FCONE FCONE);
        }
      }else base_rhs=A.transpose()*response;
      rhs_ready=true;
    }
    return base_rhs;
  }
  bool solve_sparse(const Eigen::SparseMatrix<double>& G,const Eigen::VectorXd& b,
      Eigen::VectorXd& answer,bool natural_ordering=false) {
    bool same=natural_mode==natural_ordering&&factor_outer.size()==(size_t)G.cols()+1&&factor_inner.size()==(size_t)G.nonZeros()&&
      std::equal(factor_outer.begin(),factor_outer.end(),G.outerIndexPtr())&&
      std::equal(factor_inner.begin(),factor_inner.end(),G.innerIndexPtr());
    natural_mode=natural_ordering;
    if(!same) {
      factor_outer.assign(G.outerIndexPtr(),G.outerIndexPtr()+G.cols()+1);
      factor_inner.assign(G.innerIndexPtr(),G.innerIndexPtr()+G.nonZeros());
    }
    if(natural_ordering) {
      if(!natural){natural=std::make_unique<Natural>();same=false;}
      if(!same)natural->analyzePattern(G);
      natural->factorize(G);if(natural->info()!=Eigen::Success)return false;
      answer=natural->solve(b);return natural->info()==Eigen::Success;
    }
    if(!amd){amd=std::make_unique<AMD>();same=false;}
    if(!same)amd->analyzePattern(G);
    amd->factorize(G);if(amd->info()!=Eigen::Success)return false;
    answer=amd->solve(b);return amd->info()==Eigen::Success;
  }
};
}
#endif
