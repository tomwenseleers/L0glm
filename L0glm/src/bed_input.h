#pragma once
#include <Rcpp.h>
#include <fstream>
#include <vector>
#include <cstdint>

struct BedInput {
  std::ifstream file;
  std::uint64_t stride;
  int columns;
  BedInput(const std::string& path,int n,int p) : file(path,std::ios::binary),
    stride((static_cast<std::uint64_t>(n)+3)/4),columns(p) {
    if(n<1||p<1||!file)Rcpp::stop("Invalid BED path or dimensions.");
    file.seekg(0,std::ios::end);
    if(static_cast<std::uint64_t>(file.tellg())!=3+stride*static_cast<std::uint64_t>(p))
      Rcpp::stop("BED size does not match n and p.");
    file.seekg(0);unsigned char header[3];file.read(reinterpret_cast<char*>(header),3);
    if(!file||header[0]!=0x6c||header[1]!=0x1b||header[2]!=0x01)
      Rcpp::stop("Expected variant-major PLINK 1 BED header.");
  }
  void read(int first,int count,std::vector<unsigned char>& buffer) {
    if(first<0||count<0||first>columns-count)Rcpp::stop("BED column index out of range.");
    buffer.resize(stride*count);file.clear();
    file.seekg(static_cast<std::streamoff>(3+stride*static_cast<std::uint64_t>(first)));
    file.read(reinterpret_cast<char*>(buffer.data()),static_cast<std::streamsize>(buffer.size()));
    if(!file)Rcpp::stop("Incomplete BED read.");
  }
};
