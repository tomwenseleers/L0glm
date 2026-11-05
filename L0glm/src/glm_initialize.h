#ifndef L0GLM_INITIALIZE_H
#define L0GLM_INITIALIZE_H

using DenseDesign=Eigen::Map<MatrixXd>;
using SparseDesign=Eigen::MappedSparseMatrix<double>;

inline MatrixXd sym_gram(const MatrixXd& A) {
 const int n=A.rows(),p=A.cols();MatrixXd G=MatrixXd::Zero(p,p);
 if(n&&p){double one=1,zero=0;
  F77_CALL(dsyrk)("L","T",&p,&n,&one,A.data(),&n,&zero,G.data(),&p FCONE FCONE);}
 G.triangularView<Eigen::Upper>()=G.transpose().triangularView<Eigen::Upper>();return G;
}
inline SparseMatrix<double> sym_gram(const SparseMatrix<double>& A) {return A.transpose()*A;}
inline MatrixXd right_product(const MatrixXd& A,const MatrixXd& B) {
 const int n=A.rows(),p=A.cols(),m=B.cols();MatrixXd out=MatrixXd::Zero(p,m);
 if(n&&p&&m){double one=1,zero=0;
  F77_CALL(dgemm)("T","N",&p,&m,&n,&one,A.data(),&n,B.data(),&n,&zero,out.data(),&p FCONE FCONE);}
 return out;
}
inline MatrixXd right_product(const SparseMatrix<double>& A,const MatrixXd& B) {return A.transpose()*B;}
inline MatrixXd qr_solve(const MatrixXd& A,const MatrixXd& B) {
 if(!A.cols())return MatrixXd(0,B.cols());
 auto qr=A.colPivHouseholderQr();qr.setThreshold(1e-7);return qr.solve(B);
}
struct DenseFactor {
 MatrixXd L;
 explicit DenseFactor(const MatrixXd& A):L(A) {
  const int n=L.rows();int info=0;
  if(n)F77_CALL(dpotrf)("L",&n,L.data(),&n,&info FCONE);
  if(info)throw std::runtime_error("Native dispersion pilot Cholesky failed.");
 }
 MatrixXd solve(const MatrixXd& B) const {
  MatrixXd out=B;const int n=L.rows(),m=B.cols();int info=0;
  if(n&&m)F77_CALL(dpotrs)("L",&n,&m,L.data(),&n,out.data(),&n,&info FCONE);
  if(info)throw std::runtime_error("Native dispersion pilot solve failed.");
  return out;
 }
};
struct SparseFactor {
 Eigen::SimplicialLLT<SparseMatrix<double>> factor;
 explicit SparseFactor(const SparseMatrix<double>& A) {
  factor.compute(A);if(factor.info()!=Eigen::Success)throw std::runtime_error("Native sparse dispersion pilot Cholesky failed.");
 }
 MatrixXd solve(const MatrixXd& B) const {return factor.solve(B);}
};
inline MatrixXd row_gram_native(const MatrixXd& A) {return L0glm::row_gram(A).selfadjointView<Eigen::Lower>();}
inline SparseMatrix<double> row_gram_native(const SparseMatrix<double>& A) {return A*A.transpose();}
inline void add_diagonal(MatrixXd& A,const VectorXd& d) {A.diagonal()+=d;}
inline void add_diagonal(SparseMatrix<double>& A,const VectorXd& d) {
 for(int j=0;j<A.rows();++j)A.coeffRef(j,j)+=d[j];
 A.makeCompressed();
}
inline MatrixXd append_columns(const MatrixXd& U,const MatrixXd& A) {
 MatrixXd D(A.rows(),U.cols()+A.cols());D<<U,A;return D;
}
inline SparseMatrix<double> append_columns(const MatrixXd& U,const SparseMatrix<double>& A) {
 std::vector<Eigen::Triplet<double>> entries;entries.reserve(A.nonZeros()+U.size());
 for(int j=0;j<U.cols();++j)for(int i=0;i<U.rows();++i)if(U(i,j)!=0)entries.emplace_back(i,j,U(i,j));
 for(int j=0;j<A.cols();++j)for(SparseMatrix<double>::InnerIterator it(A,j);it;++it)
  entries.emplace_back(it.row(),U.cols()+j,it.value());
 SparseMatrix<double> D(A.rows(),U.cols()+A.cols());D.setFromTriplets(entries.begin(),entries.end());return D;
}
template<class Matrix,class Factor> struct PilotOperator {
 const Matrix& A;const MatrixXd& U;Matrix D;const Matrix* design;
 std::unique_ptr<Factor> factor;MatrixXd inverse_U,schur;double tau;bool dual;
 PilotOperator(const Matrix& a,const MatrixXd& u,double t):A(a),U(u),design(&a),tau(t),dual(a.cols()>a.rows()) {
  const int n=A.rows(),p=A.cols(),q=U.cols();if(!p)return;
  if(dual) {
   Matrix K=row_gram_native(A);add_diagonal(K,VectorXd::Constant(n,tau));
   factor=std::make_unique<Factor>(K);
   if(q){inverse_U=factor->solve(U);schur=U.transpose()*inverse_U;}
  }else {
   if(q){D=append_columns(U,A);design=&D;}
   Matrix K=sym_gram(*design);VectorXd penalty=VectorXd::Constant(p+q,tau);
   if(q)penalty.head(q).setZero();
   add_diagonal(K,penalty);factor=std::make_unique<Factor>(K);
  }
 }
 MatrixXd apply(const MatrixXd& B) const {
  if(!A.cols())return B;
  if(dual) {
   MatrixXd value=factor->solve(B);
   if(U.cols())value.noalias()-=inverse_U*qr_solve(schur,U.transpose()*value);
   return tau*value;
  }
  return B-(*design)*factor->solve(right_product(*design,B));
 }
};
inline MatrixXd select_columns(const MatrixXd& A,const std::vector<int>& cols) {
 MatrixXd B(A.rows(),cols.size());for(int j=0;j<(int)cols.size();++j)B.col(j)=A.col(cols[j]);return B;
}
inline SparseMatrix<double> select_columns(const SparseMatrix<double>& A,const std::vector<int>& cols) {
 std::vector<Eigen::Triplet<double>> entries;
 for(int j=0;j<(int)cols.size();++j)for(SparseMatrix<double>::InnerIterator it(A,cols[j]);it;++it)
  entries.emplace_back(it.row(),j,it.value());
 SparseMatrix<double> B(A.rows(),cols.size());B.setFromTriplets(entries.begin(),entries.end());return B;
}
inline MatrixXd weighted_design(const DenseDesign& X,const std::vector<int>& rows,
 const std::vector<int>& cols,const VectorXd& W) {
 MatrixXd A(rows.size(),cols.size());VectorXd root(rows.size());
 for(int i=0;i<(int)rows.size();++i)root[i]=std::sqrt(W[rows[i]]);
 for(int j=0;j<(int)cols.size();++j)for(int i=0;i<(int)rows.size();++i)A(i,j)=X(rows[i],cols[j])*root[i];
 return A;
}
inline SparseMatrix<double> weighted_design(const SparseDesign& X,const std::vector<int>& rows,
 const std::vector<int>& cols,const VectorXd& W) {
 std::vector<int> map(X.rows(),-1);for(int i=0;i<(int)rows.size();++i)map[rows[i]]=i;
 std::vector<Eigen::Triplet<double>> entries;VectorXd root=W.cwiseSqrt();
 for(int j=0;j<(int)cols.size();++j)for(SparseDesign::InnerIterator it(X,cols[j]);it;++it)
  if(map[it.row()]>=0&&it.value()!=0)entries.emplace_back(map[it.row()],j,it.value()*root[it.row()]);
 SparseMatrix<double> A(rows.size(),cols.size());A.setFromTriplets(entries.begin(),entries.end());return A;
}
template<class Design> void moments(const Design& X,const std::vector<int>& rows,
 const std::vector<int>& cols,const VectorXd& W,VectorXd& ss,VectorXd& centered) {
 const int p=cols.size();ss=VectorXd::Zero(p);centered=VectorXd::Zero(p);
 double total=0;for(int i:rows)total+=W[i];
 for(int j=0;j<p;++j) {
  const double origin=rows.empty()?0:X.coeff(rows[0],cols[j]);double sum=0;
  for(int i:rows){double x=X.coeff(i,cols[j]);sum+=W[i]*(x-origin);ss[j]+=W[i]*x*x;}
  const double mean=total>0?origin+sum/total:0;
  for(int i:rows){double d=X.coeff(i,cols[j])-mean;centered[j]+=W[i]*d*d;}
 }
}
inline void moments(const SparseDesign& X,const std::vector<int>& rows,
 const std::vector<int>& cols,const VectorXd& W,VectorXd& ss,VectorXd& centered) {
 const int p=cols.size();ss=VectorXd::Zero(p);centered=VectorXd::Zero(p);
 std::vector<char> good(X.rows(),0);double total=0;for(int i:rows){good[i]=1;total+=W[i];}
 for(int j=0;j<p;++j) {
  int count=0;double first=0,sum=0,represented=0;
  for(SparseDesign::InnerIterator it(X,cols[j]);it;++it)if(good[it.row()]){if(!count)first=it.value();++count;}
  const double origin=count==(int)rows.size()?first:0;
  for(SparseDesign::InnerIterator it(X,cols[j]);it;++it)if(good[it.row()]) {
   const double w=W[it.row()],x=it.value();sum+=w*(x-origin);ss[j]+=w*x*x;represented+=w;
  }
  const double mean=total>0?origin+sum/total:0;
  for(SparseDesign::InnerIterator it(X,cols[j]);it;++it)if(good[it.row()]) {
   double d=it.value()-mean;centered[j]+=W[it.row()]*d*d;
  }
  centered[j]+=std::max(0.,total-represented)*mean*mean;
 }
}
inline void initialize_working(State& s,const Family& f,bool scalar_score) {
 const int n=s.y.size();s.W=VectorXd::Zero(n);s.z=VectorXd::Zero(n);
 for(int i=0;i<n;++i) {
  const double d=f.derivative(s.eta[i]);double w=s.w[i]*d*d/f.variance(s.mu[i]);
  double z=f.link==0?s.y[i]-s.offset[i]:s.eta[i]-s.offset[i]+(s.y[i]-s.mu[i])/d;
  if(scalar_score&&f.scalar>=0) {
   const double mu=std::max(s.mu[i],L0glmScalar::tiny);
   const auto a=L0glmScalar::derivatives(s.y[i],mu,f.scalar,f.parameter,f.theta,f.censor);
   const double information=std::max({1e-20,a.information,s.y[i]==0?1/(f.parameter*mu):0.});
   w=s.w[i]*information*d*d;z=s.eta[i]-s.offset[i]-a.gradient/(information*d);
  }
  if(s.w[i]==0)w=0;
  if(s.Wret.size()==n){s.Wret[i]=w;s.zret[i]=z;}
  if(w>0&&std::isfinite(w)&&std::isfinite(z)&&d!=0){s.W[i]=w;s.z[i]=z;}
 }
}
template<class Design> std::vector<int> eligible_columns(const Design& X,const State& s) {
 std::vector<int> cols;for(int j=0;j<X.cols();++j)if(s.allowed[j]||s.included[j])cols.push_back(j);return cols;
}
inline std::vector<int> positive_rows(const VectorXd& W) {
 std::vector<int> rows;for(int i=0;i<W.size();++i)if(W[i]>0)rows.push_back(i);return rows;
}
template<class Design> bool unit_column(const Design& X,int j,const std::vector<int>& rows) {
 if(rows.empty())return false;
 for(int i:rows)if(X.coeff(i,j)!=1)return false;
 return true;
}
template<class Design> VectorXd marginal_start(const Design& X,State& s,const Family& f) {
 const int p=X.cols();VectorXd result=VectorXd::Zero(p);
 initialize_working(s,f,true);auto rows=positive_rows(s.W);auto cols=eligible_columns(X,s);
 if(rows.empty()||cols.empty())return result;
 auto A=weighted_design(X,rows,cols,s.W);VectorXd b(rows.size());
 for(int i=0;i<(int)rows.size();++i)b[i]=s.z[rows[i]]*std::sqrt(s.W[rows[i]]);
 VectorXd info,centered;moments(X,rows,cols,s.W,info,centered);
 VectorXd rhs=right_product(A,MatrixXd(b)).col(0);
 std::vector<int> free;
 if(s.normalized) {
  for(int j=0;j<(int)cols.size();++j)if(s.included[cols[j]]&&s.lower[cols[j]]<0)free.push_back(j);
 }else {
  // Legacy row augmentation forms sums of weighted squares, not original moments.
  for(int j=0;j<(int)cols.size();++j)info[j]=A.col(j).squaredNorm();
 }
 MatrixXd U=MatrixXd(select_columns(A,free)),G=U.transpose()*U;VectorXd response=b;
 if(!free.empty()) {
  MatrixXd C=U.transpose()*A;
  response-=U*qr_solve(G,U.transpose()*b);
  info=(info.array()-(C.array()*qr_solve(G,C).array()).colwise().sum().transpose()).max(0);
  if(free.size()==1&&unit_column(X,cols[free[0]],rows))info=centered;
  rhs=right_product(A,MatrixXd(response)).col(0);
 }
 for(int j=0;j<(int)cols.size();++j) {
  int a=cols[j];double penalty=s.lambda*(s.normalized?s.scale[a]*s.scale[a]:1);
  if(s.normalized&&s.included[a])penalty=0;
  result[a]=info[j]+penalty>0?(rhs[j]+penalty*s.prior[a])/(info[j]+penalty):0;
  if(s.normalized&&info[j]==0&&s.prior[a]==0)result[a]=0;
 }
 if(!s.normalized&&s.binomial_name) {
  std::vector<int> intercept;
  for(int j=0;j<(int)cols.size();++j)if(unit_column(X,cols[j],rows))intercept.push_back(j);
  if(!intercept.empty()) {
   int k=intercept[0];double ii=info[k]+s.lambda;
   for(int j=0;j<(int)cols.size();++j)if(std::find(intercept.begin(),intercept.end(),j)==intercept.end()&&ii>0) {
    double c=0;for(int i:rows)c+=s.W[i]*X.coeff(i,cols[j]);
    double denom=info[j]+s.lambda-c*c/ii;
    result[cols[j]]=denom>0?(rhs[j]+s.lambda*s.prior[cols[j]]-c*(rhs[k]+s.lambda*s.prior[cols[k]])/ii)/denom:0;
   }
  }
 }
 for(int j:free)result[cols[j]]=0;
 if(s.gaussian_name) {
  VectorXd local(cols.size());for(int j=0;j<(int)cols.size();++j)local[j]=result[cols[j]];
  VectorXd prediction=A*local;
  if(!free.empty())prediction-=U*qr_solve(G,U.transpose()*prediction);
  const double den=prediction.squaredNorm();
  double factor=den>0?prediction.dot(response)/den:0;
  for(int j:cols)result[j]*=factor;
 }
 if(!free.empty()) {
  VectorXd local(cols.size());for(int j=0;j<(int)cols.size();++j)local[j]=result[cols[j]];
  VectorXd alpha=qr_solve(G,U.transpose()*(b-A*local)).col(0);
  for(int j=0;j<(int)free.size();++j)result[cols[free[j]]]=alpha[j];
 }
 return result;
}
inline double quantile05(std::vector<double> values) {
 if(values.empty())return 0;
 const double index=(values.size()-1)*.05;int lo=(int)index,hi=(int)std::ceil(index);
 std::nth_element(values.begin(),values.begin()+lo,values.end());double a=values[lo];
 if(hi==lo)return a;
 std::nth_element(values.begin(),values.begin()+hi,values.end());return a+(index-lo)*(values[hi]-a);
}
template<class Design> bool signal_start(const Design& X,State& s,VectorXd& result,bool require_offset=true) {
 if(!s.signal_family||s.maxit<=1||!(s.lower.array()>=0).all()||(require_offset&&s.offset.isZero()))return false;
 if(!L0glm::boundary_nonnegative(X))return false;
 auto cols=eligible_columns(X,s);std::vector<int> mandatory,selectable;
 for(int j:cols)(s.included[j]?mandatory:selectable).push_back(j);
 if(selectable.empty()||s.npositive==0)return false;
 VectorXd residual=(s.y-s.offset).cwiseMax(0);result=VectorXd::Zero(X.cols());
 std::vector<int> rows;for(int i=0;i<X.rows();++i)if(s.w[i]>0)rows.push_back(i);
 if(!mandatory.empty()) {
  MatrixXd B(X.rows(),mandatory.size());for(int j=0;j<(int)mandatory.size();++j)B.col(j)=X.col(mandatory[j]);
  VectorXd mass=B.rowwise().sum();if(!(mass.array()>0).any())return false;
  std::vector<double> vals;for(int i:rows)vals.push_back(residual[i]);double level=quantile05(vals);
  VectorXd a;bool constant=mass[0]>0&&(mass.array()-mass[0]).abs().maxCoeff()<=1e-12*std::max(1.,mass[0]);
  if(constant)a=VectorXd::Constant(mandatory.size(),level/mass[0]);
  else {
   MatrixXd weighted=B.array().colwise()*s.w.array();
   a=L0glm::nuisance_qp(B.transpose()*weighted,B.transpose()*(s.w*level),std::vector<char>(mandatory.size(),1));
  }
  for(int j=0;j<(int)mandatory.size();++j)result[mandatory[j]]=a[j];
  residual=(residual-B*a).cwiseMax(0);
 }
 double magnitude=0;for(int i:rows)magnitude+=s.w[i]*residual[i]*residual[i];magnitude=1e-8*std::sqrt(magnitude);
 for(int j:selectable) {
  double info=0,rhs=0;
  L0glm::boundary_column(X,j,[&](int i,double x){info+=s.w[i]*x*x;rhs+=s.w[i]*x*residual[i];});
  result[j]=info==0?0:std::max({rhs/std::max(DBL_MIN,info),magnitude/std::sqrt(std::max(DBL_MIN,info)),s.prior[j]});
 }
 return true;
}
template<class Design> List signal_matrix(const Design& X,NumericMatrix Y,NumericMatrix O,
 NumericMatrix W,NumericMatrix P,IntegerVector mandatory,int nthreads) {
 const int n=Y.nrow(),m=Y.ncol(),p=X.cols();MatrixXd result=MatrixXd::Zero(p,m);
 if(X.rows()!=n||O.nrow()!=n||W.nrow()!=n||P.nrow()!=p||
    !(O.ncol()==1||O.ncol()==m)||!(W.ncol()==1||W.ncol()==m)||!(P.ncol()==1||P.ncol()==m))
  stop("Invalid native signal-start dimensions.");
 std::vector<State> states(m);std::vector<char> valid(m,0);std::vector<std::string> errors(m);
 for(int k=0;k<m;++k) {
  State& s=states[k];s.y=Eigen::Map<VectorXd>(Y.begin()+n*k,n);
  s.offset=Eigen::Map<VectorXd>(O.begin()+n*(O.ncol()==1?0:k),n);
  s.w=Eigen::Map<VectorXd>(W.begin()+n*(W.ncol()==1?0:k),n);
  s.prior=Eigen::Map<VectorXd>(P.begin()+p*(P.ncol()==1?0:k),p);
  s.lower=VectorXd::Zero(p);s.included.assign(p,0);s.allowed.assign(p,1);s.signal_family=true;s.maxit=25;s.npositive=0;
  for(int j:mandatory)if(j>0&&j<=p)s.included[j-1]=1;
  for(double w:s.w)if(w>0)++s.npositive;
 }
 #ifdef _OPENMP
 #pragma omp parallel for schedule(dynamic,1) num_threads(nthreads)
 #endif
 for(int k=0;k<m;++k)try {
  VectorXd beta;valid[k]=signal_start(X,states[k],beta,false);
  if(valid[k])result.col(k)=beta;
 }catch(const std::exception& e){errors[k]=e.what();}
 for(const auto& e:errors)if(!e.empty())stop(e);
 bool informative=false;for(char v:valid)informative=informative||v;
 return List::create(_["coefficients"]=result,_["valid"]=informative);
}

template<class Matrix> double dispersion_from_design(const Matrix& A,const MatrixXd& U,const MatrixXd& B,
 bool canonical,double reference,State& s) {
 const int n=A.rows();const double target=std::min(10.,(n-U.cols())/5.);
 auto calculate=[&](double tau) {
  using Factor=std::conditional_t<std::is_same_v<Matrix,MatrixXd>,DenseFactor,SparseFactor>;
  PilotOperator<Matrix,Factor> op(A,U,tau);MatrixXd residual=op.apply(B);
  double df2=residual.rightCols(4).squaredNorm()/4,eff=df2;
  if(canonical&&(!std::isfinite(df2)||df2<target)) {
   MatrixXd repeated=op.apply(MatrixXd(residual.rightCols(4)));
   double df4=repeated.squaredNorm()/4;eff=std::isfinite(df4)&&df4>0?df2*df2/df4:0;
  }
  double phi=std::isfinite(df2)&&df2>0?residual.col(0).squaredNorm()/df2:NA_REAL;
  return std::make_pair(phi,eff);
 };
 auto value=calculate(1e-6);
 double null_scale=B.col(0).squaredNorm()/std::max(1,(int)(n-U.cols()));
 if(canonical&&(!std::isfinite(value.first)||value.first<1e-6*null_scale||value.second<target)) {
  value=calculate(.1);
 }
 double phi=value.first;if(!std::isfinite(phi)||phi<=0)phi=eps;
 return std::max(DBL_MIN,phi/reference);
}
template<class Design> double initial_dispersion(const Design& X,State& s,const Family& f,
 const std::vector<double>& probes) {
 initialize_working(s,f,false);auto rows=positive_rows(s.W);auto cols=eligible_columns(X,s);
 if(rows.empty())return 1;
 if(cols.empty())return 1;
 auto A=weighted_design(X,rows,cols,s.W);VectorXd info,centered;moments(X,rows,cols,s.W,info,centered);
 std::vector<int> free;
 if(s.normalized)for(int j=0;j<(int)cols.size();++j)if(s.included[cols[j]])free.push_back(j);
 MatrixXd U=MatrixXd(select_columns(A,free));
 if(U.cols()) {
  auto qr=U.colPivHouseholderQr();qr.setThreshold(1e-7);
  if(qr.rank()<U.cols()) {
   std::vector<int> keep;for(int j=0;j<qr.rank();++j)keep.push_back(qr.colsPermutation().indices()[j]);
   std::sort(keep.begin(),keep.end());U=select_columns(U,keep);
  }
 }
 MatrixXd G=U.transpose()*U;
 if(U.cols()) {
  MatrixXd C=U.transpose()*A;
  info=(info.array()-(C.array()*qr_solve(G,C).array()).colwise().sum().transpose()).max(0);
  if(U.cols()==1&&free.size()==1&&unit_column(X,cols[free[0]],rows))info=centered;
 }
 std::vector<int> use;
 for(int j=0;j<(int)cols.size();++j)if((!s.normalized||!s.included[cols[j]])&&std::isfinite(info[j])&&info[j]>0)use.push_back(j);
 auto scaled=select_columns(A,use);VectorXd inv(use.size());for(int j=0;j<(int)use.size();++j)inv[j]=1/std::sqrt(info[use[j]]);
 scaled=scaled*inv.asDiagonal();
 if constexpr(std::is_same_v<decltype(scaled),MatrixXd>) {
  if(U.cols())scaled-=U*qr_solve(G,U.transpose()*scaled);
 }
 const int n=rows.size();MatrixXd B(n,5);
 for(int i=0;i<n;++i){B(i,0)=s.z[rows[i]]*std::sqrt(s.W[rows[i]]);
  for(int j=0;j<4;++j)B(i,j+1)=probes[i+n*j];}
 if(U.cols())B-=U*qr_solve(G,U.transpose()*B);
 double p=s.normalized?std::max(1.,s.reference_p):(double)cols.size();
 double count=s.normalized?n:X.rows();
 double reference=std::max(1.,std::log(p)*std::log(std::log(count)));
 if(!s.normalized&&!scaled.cols())return std::max(DBL_MIN,B.col(0).squaredNorm()/n/reference);
 if constexpr(std::is_same_v<decltype(scaled),MatrixXd>) {
  if(scaled.size()>=1000000&&(scaled.array()!=0).count()<.1*scaled.size()) {
   SparseMatrix<double> sparse=scaled.sparseView();
   return dispersion_from_design(sparse,U,B,s.normalized,reference,s);
  }
 }
 return dispersion_from_design(scaled,U,B,s.normalized,reference,s);
}
template<class Design> void prepare_state(const Design& X,State& s,const Family& f,
 const std::vector<double>& probes) {
 if(s.auto_beta) {
  VectorXd beta;
  if(!signal_start(X,s,beta))beta=marginal_start(X,s,f);
  s.beta=beta.cwiseMax(s.lower).cwiseMin(s.upper);
 }
 s.beta=s.beta.cwiseMax(s.lower).cwiseMin(s.upper);
 if(s.auto_phi) {
  if(eligible_columns(X,s).empty()) {
   double value=0;for(int i=0;i<X.rows();++i)if(s.w[i]>0) {
    double mu=f.mean(s.offset[i]);value+=s.w[i]*std::pow(s.y[i]-mu,2)/f.variance(mu);
   }
   s.phi=s.npositive>0?value/s.npositive:NA_REAL;
  }else s.phi=initial_dispersion(X,s,f,probes);
  if(!std::isfinite(s.phi)||s.phi<=0)s.phi=1;
 }
 s.anchor=s.beta;s.previous_eta=predict(X,s.beta,s.offset);
 s.previous_dev=f.deviance(s.y,s.previous_eta,s.w);
 if(!std::isfinite(s.previous_dev)&&std::isfinite(f.deviance(s.y,s.offset,s.w))) {
  const VectorXd previous=s.previous_eta;
  for(int attempt=1;attempt<=60;++attempt) {
   VectorXd trial=std::ldexp(1.,-attempt)*s.beta;
   VectorXd eta=previous.allFinite()?s.offset+std::ldexp(1.,-attempt)*(previous-s.offset):predict(X,trial,s.offset);
   double value=f.deviance(s.y,eta,s.w);
   if(std::isfinite(value)){s.anchor=trial;s.previous_eta=eta;s.previous_dev=value;break;}
  }
  if(!std::isfinite(s.previous_dev)){s.anchor.setZero();s.previous_eta=s.offset;s.previous_dev=f.deviance(s.y,s.offset,s.w);}
 }
 if(s.normalized&&std::any_of(s.included.begin(),s.included.end(),[](char x){return x;}))
  s.old_prediction=predict(X,s.beta,VectorXd::Zero(X.rows()));
}
#endif
