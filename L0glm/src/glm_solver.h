#ifndef L0GLM_GLM_SOLVER_H
#define L0GLM_GLM_SOLVER_H

// This header is included inside L0glmNative, after the working-design helpers.
struct SolverControl {
  int method=0, maxit=0;
  unsigned int pivot=0, ordering=0, precond=0, svd=0, sparse_svd=1;
  double tol=0, threshold=0;
  bool rowaug=false;
  SEXP bridge=R_NilValue;
  SolverControl()=default;
  template<class T> T option(List x,const char* name,T fallback) {
    return x.containsElementNamed(name)&&!Rf_isNull(x[name])?as<T>(x[name]):fallback;
  }
  SolverControl(List x,SEXP callback):bridge(callback) {
    std::string type=option<std::string>(x,"method","chol");
    method=type=="CG"?1:type=="LSCG"?2:type=="QR"?3:type=="SVD"?4:0;
    pivot=option<unsigned int>(x,"pivot",0);ordering=option<unsigned int>(x,"ordering",0);
    precond=option<unsigned int>(x,"precond",0);svd=option<unsigned int>(x,"svd",0);
    sparse_svd=option<unsigned int>(x,"sparse.svd",1);maxit=option<int>(x,"maxit",0);
    tol=option<double>(x,"tol",0);threshold=option<double>(x,"threshold",0);
    rowaug=option<bool>(x,"rowaug",false);
  }
  bool external() const {return bridge!=R_NilValue;}
};
template<class Matrix> L0glm::SolverResult eigen_solve(const Matrix& A,const MatrixXd& b,
  const MatrixXd& start,const SolverControl& control) {
  using Ref=std::conditional_t<std::is_same_v<Matrix,MatrixXd>,L0glm::DenseRef,L0glm::SparseRef>;
  Ref design(A);L0glm::DenseRef rhs(b),guess(start);
  switch(control.method) {
    case 0:return L0glm::chol_solve(design,rhs,control.pivot,control.ordering);
    case 1:return L0glm::cg_solve(design,rhs,guess,control.tol,control.maxit,control.precond);
    case 2:return L0glm::lscg_solve(design,rhs,guess,control.tol,control.maxit,control.precond);
    case 3:return L0glm::qr_solve(design,rhs,control.pivot,control.ordering);
    case 4:return L0glm::svd_solve(design,rhs,std::is_same_v<Matrix,MatrixXd>?control.svd:control.sparse_svd,
      control.precond,control.threshold);
  }
  return L0glm::SolverResult();
}
template<class Matrix> bool eigen_coefficients(const Matrix& A,const MatrixXd& b,
  const MatrixXd& start,const SolverControl& control,State& state,MatrixXd& answer,
  bool least_squares=false) {
  auto result=eigen_solve(A,b,start,control);
  if(result.finite()) {
    if(!result.success && (control.method==1||control.method==2))++state.solver_inexact;
    answer=result.coefficients;return true;
  }
  if(control.method!=0) {
    ++state.solver_fallbacks;SolverControl fallback;
    if(least_squares) {
      Matrix gram=sym_gram(A);
      MatrixXd rhs=right_product(A,b);
      result=eigen_solve(gram,rhs,start,fallback);
    }else result=eigen_solve(A,b,start,fallback);
    if(result.finite()){answer=result.coefficients;return true;}
  }
  return false;
}
inline SparseMatrix<double> sparse_copy(const MatrixXd& A) {return A.sparseView();}
inline SparseMatrix<double> sparse_copy(const SparseMatrix<double>& A) {return A;}
template<class Matrix> SparseMatrix<double> augmented_design(const Matrix& A,const VectorXd& d) {
  SparseMatrix<double> sparse=sparse_copy(A);
  std::vector<Eigen::Triplet<double>> entries;entries.reserve(sparse.nonZeros()+d.size());
  for(int j=0;j<sparse.cols();++j)for(SparseMatrix<double>::InnerIterator it(sparse,j);it;++it)
    entries.emplace_back(it.row(),j,it.value());
  for(int j=0;j<d.size();++j)if(d[j]!=0)entries.emplace_back(A.rows()+j,j,std::sqrt(d[j]));
  SparseMatrix<double> result(A.rows()+d.size(),A.cols());
  result.setFromTriplets(entries.begin(),entries.end());return result;
}
template<class Matrix> bool selected_dual(const Matrix& A,const VectorXd& response,
  const VectorXd& penalty,VectorXd center,const SolverControl& control,State& state,VectorXd& fit) {
  const int n=A.rows(),p=A.cols();bool penalized=(penalty.array()>0).any();
  std::vector<int> free;VectorXd inv=VectorXd::Ones(p);
  for(int j=0;j<p;++j)if(penalized) {
    inv[j]=penalty[j]>0?1/std::sqrt(penalty[j]):0;
    if(penalty[j]==0){free.push_back(j);center[j]=0;}
  }else center[j]=0;
  Matrix Z=A*inv.asDiagonal();L0glm::drop_exact_zeros(Z);
  Matrix G=row_gram_native(Z);
  if(penalized)add_diagonal(G,VectorXd::Ones(n));
  MatrixXd rhs(n,free.size()+1),solution,U(n,free.size());rhs.col(0)=response-A*center;
  for(int j=0;j<(int)free.size();++j){U.col(j)=A.col(free[j]);rhs.col(j+1)=U.col(j);}
  if(!eigen_coefficients(G,rhs,MatrixXd::Zero(n,rhs.cols()),control,state,solution))return false;
  VectorXd u=solution.col(0),alpha(free.size());
  if(!free.empty()) {
    MatrixXd T=solution.rightCols(free.size()),K=U.transpose()*T;
    alpha=qr_solve(K,U.transpose()*u);u.noalias()-=T*alpha;
  }
  fit=center+inv.asDiagonal()*(Z.transpose()*u);
  for(int j=0;j<(int)free.size();++j)fit[free[j]]=alpha[j];
  return fit.allFinite();
}
template<class Matrix> bool primal_ridge(const Matrix& A,const VectorXd& response,
  const VectorXd& penalty,const VectorXd& center,const VectorXd& start,
  const SolverControl& control,State& state,VectorXd& fit) {
  Matrix G=sym_gram(A);add_diagonal(G,penalty);
  MatrixXd rhs=right_product(A,MatrixXd(response));rhs.col(0)+=penalty.cwiseProduct(center);
  MatrixXd solution;
  if(!eigen_coefficients(G,rhs,MatrixXd(start),control,state,solution))return false;
  fit=solution.col(0);return true;
}
template<class Design> VectorXd ridge_selected(const Design& X,State& state,
  const std::vector<int>& active,const SolverControl& control,int iteration) {
  auto rows=positive_rows(state.W);const int n=rows.size(),p=active.size();
  auto A=weighted_design(X,rows,active,state.W);
  VectorXd response(n),inv(p),d(p),center(p),guess(p),lower(p),upper(p);
  bool scaled=(state.scale.array()!=1).any();
  for(int i=0;i<n;++i)response[i]=state.z[rows[i]]*std::sqrt(state.W[rows[i]]);
  for(int j=0;j<p;++j) {
    const int a=active[j];inv[j]=1/state.scale[a];d[j]=state.lambdas[a]/(state.scale[a]*state.scale[a]);
    center[j]=state.prior[a]*state.scale[a];guess[j]=state.beta[a]*state.scale[a];
    lower[j]=scaled&&state.lower[a]<=-DBL_MAX?R_NegInf:state.lower[a]*state.scale[a];
    upper[j]=scaled&&state.upper[a]>=DBL_MAX?R_PosInf:state.upper[a]*state.scale[a];
  }
  A=A*inv.asDiagonal();VectorXd fit(p);bool success=false;
  if(control.external()) {
    // Dependency wrappers allocate R objects; callers execute this branch on
    // the R thread (or in independent R workers), never in OpenMP workers.
    Function solve(control.bridge);
    fit=as<VectorXd>(solve(wrap(A),wrap(response),wrap(d),wrap(center),wrap(guess),
      wrap(lower),wrap(upper),iteration));success=fit.size()==p&&fit.allFinite();
  }else if(p>n && (d.array()>=0).all() && ((d.array()==0).all() || (d.array()==0).count()<=n)) {
    success=selected_dual(A,response,d,center,control,state,fit);
  }else if(control.rowaug && (d.array()!=0).any()) {
    SparseMatrix<double> augmented=augmented_design(A,d);VectorXd rhs(n+p);rhs.head(n)=response;
    if((center.array()==0).all())rhs.tail(p).setZero();else rhs.tail(p)=d.cwiseSqrt().cwiseProduct(center);
    if(control.method>=2) {
      MatrixXd solution;success=eigen_coefficients(augmented,MatrixXd(rhs),MatrixXd(guess),control,state,solution,true);
      if(success)fit=solution.col(0);
    }else success=primal_ridge(augmented,rhs,VectorXd::Zero(p),VectorXd::Zero(p),guess,control,state,fit);
  }else if(control.rowaug && control.method>=2) {
    MatrixXd solution;success=eigen_coefficients(A,MatrixXd(response),MatrixXd(guess),control,state,solution,true);
    if(success)fit=solution.col(0);
  }else success=primal_ridge(A,response,d,center,guess,control,state,fit);
  VectorXd result=state.beta;
  if(!success){state.solver_failed=true;return result;}
  for(int j=0;j<p;++j)result[active[j]]=fit[j]/state.scale[active[j]];
  return result;
}
#endif
