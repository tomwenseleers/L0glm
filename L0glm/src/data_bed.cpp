#include "bed_input.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// Decode PLINK 1 variant-major BED: 00=2, 01=missing, 10=1, 11=0
// copies of BIM allele 1, matching PLINK's ALT convention for BED imports.
// [[Rcpp::export]]
Rcpp::NumericMatrix l0glm_bed_read_cpp(std::string path,int n,int p,
  Rcpp::IntegerVector rows,Rcpp::IntegerVector columns,int nthreads=1) {
  if(n<1||p<1||nthreads<1)Rcpp::stop("Invalid BED dimensions or thread count.");
  BedInput file(path,n,p);const std::uint64_t stride=file.stride;
  for(int row:rows)if(row<1||row>n)Rcpp::stop("BED row index out of range.");
  for(int col:columns)if(col<1||col>p)Rcpp::stop("BED column index out of range.");
  Rcpp::NumericMatrix answer(rows.size(),columns.size());
  if(!rows.size()||!columns.size())return answer;
  std::vector<unsigned char> bytes;
  const double lookup[4]={2,NA_REAL,1,0};
  const int* row=rows.begin();double* output=answer.begin();
  int workers=1;
#ifdef _OPENMP
  workers=std::min(nthreads,omp_get_num_procs());
#endif
  for(int pos=0;pos<columns.size();) {
    int count=1;
    while(pos+count<columns.size()&&columns[pos+count]==columns[pos]+count&&count<256)++count;
    file.read(columns[pos]-1,count,bytes);
    const unsigned char* packed=bytes.data();
    const int nr=rows.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(workers) if(workers>1 && count>1)
#endif
    for(int j=0;j<count;++j) {
      const unsigned char* source=packed+j*stride;
      double* destination=output+static_cast<std::size_t>(pos+j)*nr;
      for(int i=0;i<nr;++i) {
        const int r=row[i]-1;
        destination[i]=lookup[(source[r/4]>>(2*(r%4)))&3];
      }
    }
    pos+=count;Rcpp::checkUserInterrupt();
  }
  return answer;
}
