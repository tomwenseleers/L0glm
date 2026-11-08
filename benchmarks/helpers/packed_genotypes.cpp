#include <Rcpp.h>
#include <fstream>
#include <vector>
#include <cstdint>

// Synthetic independent genotypes for storage/throughput experiments.
// BIM allele-1 counts: P(0)=1/4, P(1)=1/4, P(2)=1/2; no missing calls.
// [[Rcpp::export]]
void benchmarkBedData(std::string path,int n,int p,std::uint32_t seed=101) {
  if(n<1||p<1||!seed)Rcpp::stop("Invalid generator parameters.");
  if(std::ifstream(path).good())Rcpp::stop("The output file already exists.");
  std::ofstream file(path,std::ios::binary);
  if(!file)Rcpp::stop("Cannot create BED file.");
  const unsigned char magic[3]={0x6c,0x1b,0x01};
  file.write(reinterpret_cast<const char*>(magic),3);
  const std::size_t stride=(n+3ULL)/4;
  std::vector<unsigned char> buffer(stride*128);
  for(int first=0;first<p;first+=128) {
    const int count=std::min(128,p-first);
    for(std::size_t k=0;k<stride*count;++k) {
      seed^=seed<<13;seed^=seed>>17;seed^=seed<<5;
      unsigned int b=seed&255;
      b&=~((b&0x55)&~((b>>1)&0x55));
      buffer[k]=static_cast<unsigned char>(b);
    }
    if(n%4)for(int j=0;j<count;++j)buffer[stride*(j+1)-1]&=(1u<<(2*(n%4)))-1;
    file.write(reinterpret_cast<const char*>(buffer.data()),stride*count);
    if(first%4096==0)Rcpp::checkUserInterrupt();
  }
  if(!file)Rcpp::stop("BED write failed.");
}
