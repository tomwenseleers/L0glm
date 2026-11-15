#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <RcppEigen.h>
#ifdef _WIN32
#include <windows.h>
#elif defined(__linux__)
#include <unistd.h>
#elif defined(__APPLE__)
#include <mach/mach.h>
#include <sys/sysctl.h>
#endif
#include "memory_budget.h"
namespace L0glm {
double available_memory() {
#ifdef _WIN32
  MEMORYSTATUSEX status;status.dwLength=sizeof(status);
  if(GlobalMemoryStatusEx(&status))return (double)status.ullAvailPhys;
#elif defined(__linux__)
  long pages=sysconf(_SC_AVPHYS_PAGES),size=sysconf(_SC_PAGESIZE);
  if(pages>0&&size>0)return (double)pages*size;
#elif defined(__APPLE__)
  mach_msg_type_number_t count=HOST_VM_INFO64_COUNT;vm_statistics64_data_t vm;
  vm_size_t size=0;
  if(host_page_size(mach_host_self(),&size)==KERN_SUCCESS&&
    host_statistics64(mach_host_self(),HOST_VM_INFO64,(host_info64_t)&vm,&count)==KERN_SUCCESS)
    return (double)(vm.free_count+vm.inactive_count)*size;
#endif
  return 4.*1024*1024*1024;
}
}
// [[Rcpp::export(rng = false)]]
double solver_workspace_budget_cpp(double requested=-1,int nthreads=1) {
  return L0glm::workspace_budget(requested,nthreads);
}
// [[Rcpp::export(rng = false)]]
double solver_workspace_estimate_cpp(SEXP design,bool gram=true,bool rows=false) {
  if(Rf_isS4(design)) {
    auto X=Rcpp::as<Eigen::MappedSparseMatrix<double>>(design);
    return gram?L0glm::gram_workspace(X,rows):L0glm::matrix_workspace(X);
  }
  Rcpp::NumericMatrix mat(design);Eigen::Map<Eigen::MatrixXd> X(mat.begin(),mat.nrow(),mat.ncol());
  return gram?L0glm::gram_workspace(X,rows):L0glm::matrix_workspace(X);
}
// [[Rcpp::export(rng = false)]]
Rcpp::List ridge_matrixfree_cpp(SEXP design,Eigen::VectorXd y,Eigen::VectorXd penalty,
  Eigen::VectorXd center,Eigen::VectorXd start,double tol=1e-10,int maxit=0,bool dual=false) {
  auto validate=[&](int n,int p) {if(y.size()!=n||penalty.size()!=p||center.size()!=p||start.size()!=p||
     !y.allFinite()||!penalty.allFinite()||(penalty.array()<0).any()||
     !center.allFinite()||!start.allFinite())Rcpp::stop("Invalid matrix-free ridge inputs.");};
  L0glm::MatrixFreeResult result;
  if(Rf_isS4(design)) {
    auto X=Rcpp::as<Eigen::MappedSparseMatrix<double>>(design);
    validate(X.rows(),X.cols());
    result=dual?L0glm::matrixfree_dual(X,y,penalty,center,tol,maxit):
      L0glm::matrixfree_ridge(X,y,penalty,center,start,tol,maxit);
  }else {
    Rcpp::NumericMatrix mat(design);Eigen::Map<Eigen::MatrixXd> X(mat.begin(),mat.nrow(),mat.ncol());
    validate(X.rows(),X.cols());
    result=dual?L0glm::matrixfree_dual(X,y,penalty,center,tol,maxit):
      L0glm::matrixfree_ridge(X,y,penalty,center,start,tol,maxit);
  }
  return Rcpp::List::create(Rcpp::_["coefficients"]=result.coefficients,Rcpp::_["success"]=result.success,
    Rcpp::_["iterations"]=result.iterations,Rcpp::_["error"]=result.error);
}
