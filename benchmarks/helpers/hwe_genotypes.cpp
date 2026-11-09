#include <Rcpp.h>
#include <fstream>
#include <vector>
#include <cstdint>
#include <algorithm>

// Independent diploid variants in Hardy-Weinberg equilibrium. The simulated
// minor allele is BIM allele 1, matching l0glmData's BED coefficient convention.
// [[Rcpp::export]]
Rcpp::NumericVector benchmarkHweBed(std::string path,int n,int p,
                                  std::uint32_t seed=101) {
  if(n<1||p<1||!seed)Rcpp::stop("Invalid generator dimensions or seed.");
  if(std::ifstream(path).good())Rcpp::stop("The output file already exists.");
  std::ofstream file(path,std::ios::binary);
  if(!file)Rcpp::stop("Cannot create BED file.");
  const unsigned char magic[3]={0x6c,0x1b,0x01};
  file.write(reinterpret_cast<const char*>(magic),3);
  const std::size_t stride=(n+3ULL)/4;
  std::vector<unsigned char> bytes(stride*128);
  Rcpp::NumericVector maf(p);
  auto uniform=[&]() {
    seed^=seed<<13;seed^=seed>>17;seed^=seed<<5;
    return (static_cast<double>(seed)+.5)/4294967296.;
  };
  for(int first=0;first<p;first+=128) {
    const int count=std::min(128,p-first);
    std::fill(bytes.begin(),bytes.begin()+stride*count,0);
    for(int j=0;j<count;++j) {
      const double frequency=.1+.4*uniform();maf[first+j]=frequency;
      for(int i=0;i<n;++i) {
        const int dosage=(uniform()<frequency)+(uniform()<frequency);
        const unsigned int code=dosage==0?3:(dosage==1?2:0);
        bytes[j*stride+i/4]|=code<<(2*(i%4));
      }
    }
    file.write(reinterpret_cast<const char*>(bytes.data()),stride*count);
    if(first%4096==0)Rcpp::checkUserInterrupt();
  }
  if(!file)Rcpp::stop("BED write failed.");
  return maf;
}
