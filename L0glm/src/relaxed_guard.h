#ifndef L0GLM_RELAXED_GUARD_H
#define L0GLM_RELAXED_GUARD_H
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <limits>
namespace L0glm {
inline bool relaxed_full_rank(Eigen::SparseMatrix<double> A) {
 const int p=A.cols();if(!p)return true;
 Eigen::VectorXd norm=Eigen::VectorXd::Zero(p);
 for(int j=0;j<p;++j)for(Eigen::SparseMatrix<double>::InnerIterator it(A,j);it;++it)norm[j]+=it.value()*it.value();
 for(int j=0;j<p;++j) {
  if(!(norm[j]>0)||!std::isfinite(norm[j]))return false;
  const double inverse=1/std::sqrt(norm[j]);
  for(Eigen::SparseMatrix<double>::InnerIterator it(A,j);it;++it)it.valueRef()*=inverse;
 }
 Eigen::SparseMatrix<double> G=A.transpose()*A;
 Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> factor;factor.compute(G);
 if(factor.info()!=Eigen::Success)return false;
 Eigen::SparseMatrix<double> L=factor.matrixL();
 double minimum=1,maximum=0;for(int j=0;j<p;++j){minimum=std::min(minimum,std::abs(L.coeff(j,j)));maximum=std::max(maximum,std::abs(L.coeff(j,j)));}
 return minimum>std::sqrt(p*std::numeric_limits<double>::epsilon())*maximum;
}
inline bool relaxed_full_rank(Eigen::MatrixXd A) {
 const int p=A.cols();if(!p)return true;
 if(A.size()>=1000000&&(A.array()!=0).count()<.1*A.size())return relaxed_full_rank(Eigen::SparseMatrix<double>(A.sparseView()));
 for(int j=0;j<p;++j){double norm=A.col(j).norm();if(!(norm>0)||!std::isfinite(norm))return false;A.col(j)/=norm;}
 Eigen::MatrixXd G=A.transpose()*A;Eigen::LLT<Eigen::MatrixXd> factor(G);
 if(factor.info()!=Eigen::Success)return false;
 Eigen::VectorXd diagonal=factor.matrixL().toDenseMatrix().diagonal();
 return diagonal.minCoeff()>std::sqrt(p*std::numeric_limits<double>::epsilon())*diagonal.maxCoeff();
}
inline bool relaxed_identifiable(const Eigen::Map<Eigen::MatrixXd>& X,
 const std::vector<int>& rows,const std::vector<int>& columns) {
 if(columns.size()>rows.size())return false;
 Eigen::MatrixXd A(rows.size(),columns.size());
 for(int j=0;j<(int)columns.size();++j)for(int i=0;i<(int)rows.size();++i)A(i,j)=X(rows[i],columns[j]);
 return relaxed_full_rank(std::move(A));
}
inline bool relaxed_identifiable(const Eigen::MappedSparseMatrix<double>& X,
 const std::vector<int>& rows,const std::vector<int>& columns) {
 if(columns.size()>rows.size())return false;
 std::vector<int> index(X.rows(),-1);for(int i=0;i<(int)rows.size();++i)index[rows[i]]=i;
 std::vector<Eigen::Triplet<double>> entries;
 for(int j=0;j<(int)columns.size();++j)for(Eigen::MappedSparseMatrix<double>::InnerIterator it(X,columns[j]);it;++it)
  if(index[it.row()]>=0&&it.value()!=0)entries.emplace_back(index[it.row()],j,it.value());
 Eigen::SparseMatrix<double> A(rows.size(),columns.size());A.setFromTriplets(entries.begin(),entries.end());
 return relaxed_full_rank(std::move(A));
}
// Positive weights change information magnitude, but not structural rank.
// Reuse the check when channels have the same eligible terms and observation mask.
struct RelaxedDesignGuard {
 std::vector<int> last_rows,last_columns;bool ready=false,valid=false;
 template<class Design> bool check(const Design& X,const std::vector<int>& rows,const std::vector<int>& columns) {
  if(!ready||rows!=last_rows||columns!=last_columns){valid=relaxed_identifiable(X,rows,columns);last_rows=rows;last_columns=columns;ready=true;}
  return valid;
 }
};
}
#endif
