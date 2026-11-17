#ifndef L0GLM_SPARSE_WORKSPACE_H
#define L0GLM_SPARSE_WORKSPACE_H
#include <algorithm>
#include <limits>
#include <vector>

namespace L0glm {
// Elimination-tree reach counts LLT entries without allocating the numerical
// factor. This is the symbolic graph used by Eigen's SimplicialCholesky.
inline double symbolic_factor_entries(const Eigen::SparseMatrix<double>& upper,double limit) {
  const int n=upper.cols();std::vector<int> parent(n,-1),seen(n,-1);
  double count=n;
  for(int k=0;k<n;++k) {
    seen[k]=k;
    for(Eigen::SparseMatrix<double>::InnerIterator it(upper,k);it;++it)if(it.row()<k) {
      int i=it.row();
      while(seen[i]!=k) {
        if(parent[i]<0)parent[i]=k;
        seen[i]=k;++count;if(count>limit)return count;
        i=parent[i];
      }
    }
  }
  return count;
}

// Refine a failed coarse guard using unique structural Gram entries and
// symbolic fill. Bound the probe itself before any potentially large factor
// allocation; count exact zeroes as absent, as the numerical solver does.
template<class Matrix> double sparse_gram_workspace(const Matrix& input,bool rows,double budget) {
  using Sp=Eigen::SparseMatrix<double>;
  const double rejected=std::numeric_limits<double>::infinity();
  if(!(budget>0))return rejected;
  const int p=rows?input.rows():input.cols();
  const double design=24.*input.nonZeros()+48.*(p+1.);
  if(design>budget)return rejected;
  Sp X;if(rows)X=input.transpose();else X=input;
  X.prune(0.);X.makeCompressed();
  Eigen::SparseMatrix<double,Eigen::RowMajor> byrow=X;
  std::vector<int> mark(p,-1),indices;
  std::vector<Eigen::Triplet<double>> entries;
  const double cap=(budget-design)/144.;
  if(p>cap)return rejected;
  for(int j=0;j<p;++j) {
    indices.clear();mark[j]=j;indices.push_back(j);
    for(Sp::InnerIterator it(X,j);it;++it)
      for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator jt(byrow,it.row());jt;++jt) {
        const int k=jt.col();
        if(k<=j&&mark[k]!=j){mark[k]=j;indices.push_back(k);}
      }
    if(entries.size()+indices.size()>cap)return rejected;
    std::sort(indices.begin(),indices.end());
    for(int i:indices)entries.emplace_back(i,j,1.);
  }
  Sp upper(p,p);upper.setFromTriplets(entries.begin(),entries.end());
  const double gram=2.*upper.nonZeros()-p;
  const double fixed=48*gram+48*(p+1.)+24*X.nonZeros();
  if(fixed>budget)return rejected;
  const double fill_limit=(budget-fixed)/24.;
  const double natural=symbolic_factor_entries(upper,fill_limit);
  if(natural>fill_limit)return rejected;
  Sp full=upper.selfadjointView<Eigen::Upper>();
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic,int> inverse;
  Eigen::AMDOrdering<int> ordering;ordering(full,inverse);
  Sp reordered(p,p);
  reordered.selfadjointView<Eigen::Upper>()=
    upper.selfadjointView<Eigen::Upper>().twistedBy(inverse.inverse());
  const double amd=symbolic_factor_entries(reordered,fill_limit);
  if(amd>fill_limit)return rejected;
  // Cover either supported ordering, Gram copies, the numeric factor and
  // work vectors. This is deliberately larger than their stored object sizes.
  return fixed+24*std::max(natural,amd);
}
}
#endif
