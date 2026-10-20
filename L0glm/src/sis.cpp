// Efficient-score screening without design-sized squares or channel Hessians.
#include <RcppEigen.h>
#include <queue>
#include <algorithm>
#include <exception>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {
using ConstMatrix=Eigen::Map<const MatrixXd>;
using ConstVector=Eigen::Map<const VectorXd>;
ConstMatrix sis_matrix(SEXP x) {
  NumericMatrix value(x);
  return ConstMatrix(value.begin(),value.nrow(),value.ncol());
}
ConstVector sis_vector(SEXP x) {
  NumericVector value(x);
  return ConstVector(value.begin(),value.size());
}

// R objects are unpacked on the main thread; workers use read-only Eigen views.
struct SisContext {
  ConstMatrix score, weight, nuisance, correction, probability, means, basis;
  ConstVector inverse, e, prefix;
  const int *rows, *groups;
  int n, r, q;
  bool cox, multinomial, one_sided;
  SisContext(const List& x) : score(sis_matrix(x["score"])), weight(sis_matrix(x["weight"])),
    nuisance(sis_matrix(x["nuisance"])), correction(sis_matrix(x["correction"])),
    probability(sis_matrix(x["probability"])), means(sis_matrix(x["means"])), basis(sis_matrix(x["basis"])),
    inverse(sis_vector(x["inverse"])), e(sis_vector(x["e"])), prefix(sis_vector(x["prefix"])),
    rows(INTEGER(x["rows"])), groups(INTEGER(x["groups"])), n(score.rows()), r(score.cols()),
    q(nuisance.cols()), cox(as<bool>(x["cox"])),
    multinomial(as<bool>(x["multinomial"])), one_sided(as<bool>(x["one_sided"])) {}
};

struct SisFeature {
  const SisContext& c;
  VectorXd U, diagonal, classes;
  MatrixXd I, A, class_cross;
  double risk=0, risk_correction=0;
  int last_group=0;
  SisFeature(const SisContext& context) : c(context), U(VectorXd::Zero(c.r)),
    diagonal(VectorXd::Zero(c.r)), classes(VectorXd::Zero(c.probability.cols())),
    I(MatrixXd::Zero(c.multinomial?c.r:0,c.multinomial?c.r:0)),
    A(MatrixXd::Zero(c.r,c.multinomial?c.q*c.r:c.q)),
    class_cross(MatrixXd::Zero(c.probability.cols(),c.multinomial?c.q:0)) {}
  void reset() {
    U.setZero();diagonal.setZero();classes.setZero();I.setZero();A.setZero();class_cross.setZero();
    risk=0;risk_correction=0;last_group=0;
  }
  void add(int i, double x) {
    if (!std::isfinite(x)) throw std::runtime_error("X must contain only finite values.");
    if (x==0) return;
    for (int a=0;a<c.r;++a) U[a]+=x*c.score(i,a);
    if (c.multinomial) {
      const double h=x*c.weight(i,0), h2=x*h;
      for (int k=0;k<classes.size();++k) {
        classes[k]+=h2*c.probability(i,k);
        for (int j=0;j<c.q;++j) class_cross(k,j)+=h*c.nuisance(i,j)*c.probability(i,k);
      }
      for (int a=0;a<c.r;++a) for (int b=0;b<c.r;++b) {
        const double product=c.means(i,a)*c.means(i,b);
        I(a,b)-=h2*product;
        for (int j=0;j<c.q;++j) A(a,j+b*c.q)-=h*c.nuisance(i,j)*product;
      }
    } else {
      for (int a=0;a<c.r;++a) {
        diagonal[a]+=x*x*c.weight(i,a);
        for (int j=0;j<c.q;++j) A(a,j)+=x*c.nuisance(i,j)*(c.cox?1:c.weight(i,a));
      }
      if (c.cox) {
        int group=c.groups[i]-1;
        risk_correction+=risk*risk*(c.prefix[group]-c.prefix[last_group]);
        risk+=c.e[i]*x; last_group=group;
      }
    }
  }
  std::pair<double,int> finish() {
    double statistic=0;
    if (c.multinomial) {
      ConstMatrix Q(c.basis.data(),c.basis.rows(),c.r);
      I+=Q.transpose()*classes.asDiagonal()*Q;
      const double cutoff=I.norm()*1e-12;
      for (int j=0;j<c.q;++j) {
        MatrixXd block=Q.transpose()*class_cross.col(j).asDiagonal()*Q;
        for (int a=0;a<c.r;++a) for (int b=0;b<c.r;++b) A(a,j+b*c.q)+=block(a,b);
      }
      if (c.q) {
        ConstMatrix G(c.inverse.data(),c.q*c.r,c.q*c.r);
        ConstVector correction(c.correction.data(),c.q*c.r);
        U-=A*correction; I-=A*G*A.transpose();
      }
      Eigen::SelfAdjointEigenSolver<MatrixXd> eig((I+I.transpose())/2);
      if (eig.info()!=Eigen::Success) throw std::runtime_error("Screening information eigensolve failed.");
      VectorXd projection=eig.eigenvectors().transpose()*U;
      for (int j=0;j<c.r;++j) if (eig.eigenvalues()[j]>cutoff)
        statistic+=projection[j]*projection[j]/eig.eigenvalues()[j];
    } else {
      const VectorXd raw_information=diagonal;
      if (c.cox) diagonal[0]-=risk_correction+risk*risk*(c.prefix[c.prefix.size()-1]-c.prefix[last_group]);
      for (int a=0;a<c.r;++a) {
        const double cutoff=std::abs(raw_information[a])*1e-12;
        if (c.q) {
          ConstMatrix G(c.inverse.data()+a*c.q*c.q,c.q,c.q);
          ConstVector correction(c.correction.data()+a*c.q,c.q);
          U[a]-=A.row(a).dot(correction);
          diagonal[a]-=(A.row(a)*G*A.row(a).transpose())(0,0);
        }
        const double u=c.one_sided?std::max(0.0,U[a]):U[a];
        if (diagonal[a]>cutoff) statistic+=u*u/diagonal[a];
      }
    }
    return {statistic,c.r==1?(U[0]>0)-(U[0]<0):NA_INTEGER};
  }
};

// The common single-response/null-intercept case needs only scalar reductions.
std::pair<double,int> sis_scalar(const SisContext& c,const double* x) {
  double U=0,I=0,A=0;
  for (int i=0;i<c.n;++i) {
    const double v=x[c.rows[i]-1];
    if (!std::isfinite(v)) throw std::runtime_error("X must contain only finite values.");
    U+=v*c.score(i,0); I+=v*v*c.weight(i,0);
    if (c.q) A+=v*c.weight(i,0)*c.nuisance(i,0);
  }
  const double cutoff=I*1e-12;
  if (c.q) { U-=A*c.correction(0,0); I-=A*A*c.inverse[0]; }
  const double u=c.one_sided?std::max(0.0,U):U;
  return {I>cutoff?u*u/I:0,(U>0)-(U<0)};
}

struct SisResults {
  using Entry=std::pair<double,int>;
  using Heap=std::priority_queue<Entry,std::vector<Entry>,std::greater<Entry>>;
  std::vector<Heap> heaps;
  NumericVector statistics;
  IntegerVector signs;
  double* statistic_data;
  int* sign_data;
  int count;
  bool full;
  SisResults(int p,int k,bool keep,int threads) : heaps(threads),statistics(keep?p:0),signs(keep?p:0),
    statistic_data(statistics.begin()),sign_data(signs.begin()),count(std::min(p,k)),full(keep) {}
  void retain(Heap& heap,Entry entry) {
    if (heap.size()<static_cast<std::size_t>(count)) heap.push(entry);
    else if (entry>heap.top()) {heap.pop();heap.push(entry);}
  }
  void add(int thread,int j,const std::pair<double,int>& result) {
    if (!std::isfinite(result.first)||result.first<0) throw std::runtime_error("Undefined screening statistic.");
    if (full) { statistic_data[j]=result.first; sign_data[j]=result.second; }
    retain(heaps[thread],{result.first,-j});
  }
  List finish(int threads) {
    Heap heap;
    for (auto& local:heaps) while (!local.empty()) {retain(heap,local.top());local.pop();}
    IntegerVector indices(count);NumericVector top(count);
    for (int i=count-1;i>=0;--i) {indices[i]=1-heap.top().second;top[i]=heap.top().first;heap.pop();}
    return List::create(_["indices"]=indices,_["top_statistics"]=top,
      _["nthreads"]=threads,_["statistics"]=full?static_cast<SEXP>(statistics):R_NilValue,
      _["sign"]=full?static_cast<SEXP>(signs):R_NilValue);
  }
};
}

// [[Rcpp::export]]
List sis_score_cpp(SEXP design,List context,int nsis,bool full,int nthreads=1) {
  SisContext c(context);
  if (nsis<1) stop("nsis must be positive.");
  if (nthreads<1) stop("nthreads must be positive.");
  const bool sparse=Rf_isS4(design);
  const double* values;
  const int *ptr=nullptr,*row=nullptr;
  int p,n;
  if (sparse) {
    S4 X(design);IntegerVector dims=X.slot("Dim");
    n=dims[0];p=dims[1];values=REAL(X.slot("x"));
  } else {
    NumericMatrix X(design);n=X.nrow();p=X.ncol();values=X.begin();
  }
  if (n!=c.n || p<1) stop("Invalid screening design dimensions.");
  if (sparse) {S4 X(design);ptr=INTEGER(X.slot("p"));row=INTEGER(X.slot("i"));}
  int threads=1;
#ifdef _OPENMP
  threads=std::min(nthreads,std::min(p,omp_get_num_procs()));
#endif
  SisResults result(p,nsis,full,threads);
  std::vector<int> inverse_rows(c.n);
  for (int i=0;i<c.n;++i) inverse_rows[c.rows[i]-1]=i;
  std::vector<std::exception_ptr> errors(threads);
  const bool scalar=!c.cox&&!c.multinomial&&c.r==1&&c.q<=1;
  int used=1;
  // Bounded batches allow main-thread interrupt/error checks between regions.
  const int batch=1024*std::min(threads,64);
  for (int first=0;first<p;) {
    const int last=first+std::min(batch,p-first);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads) if(threads>1)
#endif
    {
      int thread=0,team=1;
#ifdef _OPENMP
      thread=omp_get_thread_num();team=omp_get_num_threads();
#endif
      if (thread==0) used=std::max(used,team);
      const int begin=first+static_cast<long long>(last-first)*thread/team;
      const int end=first+static_cast<long long>(last-first)*(thread+1)/team;
      try {
        SisFeature feature(c);
        std::vector<std::pair<int,double>> ordered;
        for (int j=begin;j<end;++j) {
          if (!sparse && scalar) result.add(thread,j,sis_scalar(c,values+static_cast<R_xlen_t>(j)*n));
          else {
            feature.reset();
            if (sparse && c.cox) {
              ordered.clear();ordered.reserve(ptr[j+1]-ptr[j]);
              for (int k=ptr[j];k<ptr[j+1];++k) ordered.emplace_back(inverse_rows[row[k]],values[k]);
              std::sort(ordered.begin(),ordered.end());
              for (const auto& value:ordered) feature.add(value.first,value.second);
            } else if (sparse) {
              for (int k=ptr[j];k<ptr[j+1];++k) feature.add(inverse_rows[row[k]],values[k]);
            } else for (int i=0;i<c.n;++i) feature.add(i,values[static_cast<R_xlen_t>(j)*n+c.rows[i]-1]);
            result.add(thread,j,feature.finish());
          }
        }
      } catch (...) {errors[thread]=std::current_exception();}
    }
    for (const auto& error:errors) if (error) std::rethrow_exception(error);
    checkUserInterrupt();first=last;
  }
  return result.finish(used);
}
