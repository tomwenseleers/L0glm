// Efficient-score screening without design-sized squares or channel Hessians.
#include <RcppEigen.h>
#include <queue>
#include <algorithm>
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {
struct SisContext {
  NumericMatrix score, weight, nuisance, correction, probability, means, basis;
  NumericVector inverse, e, prefix;
  IntegerVector rows, groups;
  int n, r, q;
  bool cox, multinomial, one_sided;
  SisContext(const List& x) : score(x["score"]), weight(x["weight"]),
    nuisance(x["nuisance"]), correction(x["correction"]),
    probability(x["probability"]), means(x["means"]), basis(x["basis"]),
    inverse(x["inverse"]), e(x["e"]), prefix(x["prefix"]),
    rows(x["rows"]), groups(x["groups"]), n(score.nrow()), r(score.ncol()),
    q(nuisance.ncol()), cox(as<bool>(x["cox"])),
    multinomial(as<bool>(x["multinomial"])), one_sided(as<bool>(x["one_sided"])) {}
};

struct SisFeature {
  SisContext& c;
  VectorXd U, diagonal, classes;
  MatrixXd I, A, class_cross;
  double risk=0, risk_correction=0;
  int last_group=0;
  SisFeature(SisContext& context) : c(context), U(VectorXd::Zero(c.r)),
    diagonal(VectorXd::Zero(c.r)), classes(VectorXd::Zero(c.probability.ncol())),
    I(MatrixXd::Zero(c.multinomial?c.r:0,c.multinomial?c.r:0)),
    A(MatrixXd::Zero(c.r,c.multinomial?c.q*c.r:c.q)),
    class_cross(MatrixXd::Zero(c.probability.ncol(),c.multinomial?c.q:0)) {}
  void add(int i, double x) {
    if (!R_finite(x)) stop("X must contain only finite values.");
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
      Eigen::Map<MatrixXd> Q(c.basis.begin(),c.basis.nrow(),c.r);
      I+=Q.transpose()*classes.asDiagonal()*Q;
      const double cutoff=I.norm()*1e-12;
      for (int j=0;j<c.q;++j) {
        MatrixXd block=Q.transpose()*class_cross.col(j).asDiagonal()*Q;
        for (int a=0;a<c.r;++a) for (int b=0;b<c.r;++b) A(a,j+b*c.q)+=block(a,b);
      }
      if (c.q) {
        Eigen::Map<MatrixXd> G(c.inverse.begin(),c.q*c.r,c.q*c.r);
        Eigen::Map<VectorXd> correction(c.correction.begin(),c.q*c.r);
        U-=A*correction; I-=A*G*A.transpose();
      }
      Eigen::SelfAdjointEigenSolver<MatrixXd> eig((I+I.transpose())/2);
      if (eig.info()!=Eigen::Success) stop("Screening information eigensolve failed.");
      VectorXd projection=eig.eigenvectors().transpose()*U;
      for (int j=0;j<c.r;++j) if (eig.eigenvalues()[j]>cutoff)
        statistic+=projection[j]*projection[j]/eig.eigenvalues()[j];
    } else {
      const VectorXd raw_information=diagonal;
      if (c.cox) diagonal[0]-=risk_correction+risk*risk*(c.prefix[c.prefix.size()-1]-c.prefix[last_group]);
      for (int a=0;a<c.r;++a) {
        const double cutoff=std::abs(raw_information[a])*1e-12;
        if (c.q) {
          Eigen::Map<MatrixXd> G(c.inverse.begin()+a*c.q*c.q,c.q,c.q);
          Eigen::Map<VectorXd> correction(c.correction.begin()+a*c.q,c.q);
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
std::pair<double,int> sis_scalar(SisContext& c,const double* x) {
  double U=0,I=0,A=0;
  for (int i=0;i<c.n;++i) {
    const double v=x[c.rows[i]-1];
    if (!R_finite(v)) stop("X must contain only finite values.");
    U+=v*c.score[i]; I+=v*v*c.weight[i];
    if (c.q) A+=v*c.weight[i]*c.nuisance[i];
  }
  const double cutoff=I*1e-12;
  if (c.q) { U-=A*c.correction[0]; I-=A*A*c.inverse[0]; }
  const double u=c.one_sided?std::max(0.0,U):U;
  return {I>cutoff?u*u/I:0,(U>0)-(U<0)};
}

struct SisResults {
  using Entry=std::pair<double,int>;
  std::priority_queue<Entry,std::vector<Entry>,std::greater<Entry>> heap;
  NumericVector statistics;
  IntegerVector signs;
  int count;
  bool full;
  SisResults(int p,int k,bool keep) : statistics(keep?p:0), signs(keep?p:0),count(std::min(p,k)),full(keep) {}
  void add(int j,const std::pair<double,int>& result) {
    if (!std::isfinite(result.first)||result.first<0) stop("Undefined screening statistic.");
    if (full) { statistics[j]=result.first; signs[j]=result.second; }
    Entry entry={result.first,-j};
    if (heap.size()<static_cast<std::size_t>(count)) heap.push(entry);
    else if (entry>heap.top()) {heap.pop();heap.push(entry);}
  }
  List finish() {
    IntegerVector indices(count);NumericVector top(count);
    for (int i=count-1;i>=0;--i) {indices[i]=1-heap.top().second;top[i]=heap.top().first;heap.pop();}
    return List::create(_["indices"]=indices,_["top_statistics"]=top,
      _["statistics"]=full?static_cast<SEXP>(statistics):R_NilValue,
      _["sign"]=full?static_cast<SEXP>(signs):R_NilValue);
  }
};
}

// [[Rcpp::export]]
List sis_score_cpp(SEXP design,List context,int nsis,bool full) {
  SisContext c(context);
  if (nsis<1) stop("nsis must be positive.");
  if (Rf_isS4(design)) {
    S4 X(design);IntegerVector dims=X.slot("Dim"),ptr=X.slot("p"),row=X.slot("i");
    NumericVector values=X.slot("x");
    SisResults result(dims[1],nsis,full);
    IntegerVector inverse_rows(c.n);for (int i=0;i<c.n;++i) inverse_rows[c.rows[i]-1]=i;
    std::vector<std::pair<int,double>> ordered;
    for (int j=0;j<dims[1];++j) {
      SisFeature feature(c);
      if (c.cox) {
        ordered.clear();ordered.reserve(ptr[j+1]-ptr[j]);
        for (int k=ptr[j];k<ptr[j+1];++k) ordered.emplace_back(inverse_rows[row[k]],values[k]);
        std::sort(ordered.begin(),ordered.end());
        for (const auto& value:ordered) feature.add(value.first,value.second);
      } else for (int k=ptr[j];k<ptr[j+1];++k) feature.add(inverse_rows[row[k]],values[k]);
      result.add(j,feature.finish());if (j%1024==0) checkUserInterrupt();
    }
    return result.finish();
  }
  NumericMatrix X(design);SisResults result(X.ncol(),nsis,full);
  for (int j=0;j<X.ncol();++j) {
    if (!c.cox&&!c.multinomial&&c.r==1&&c.q<=1)
      result.add(j,sis_scalar(c,X.begin()+static_cast<R_xlen_t>(j)*c.n));
    else {
      SisFeature feature(c);
      for (int i=0;i<c.n;++i) feature.add(i,X(c.rows[i]-1,j));
      result.add(j,feature.finish());
    }
    if (j%1024==0) checkUserInterrupt();
  }
  return result.finish();
}
