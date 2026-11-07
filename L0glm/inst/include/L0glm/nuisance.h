#ifndef L0GLM_NUISANCE_H
#define L0GLM_NUISANCE_H
#include <vector>
#include <limits>
#include <cmath>
namespace L0glm {
inline Eigen::VectorXd nuisance_qp(const Eigen::MatrixXd& G,const Eigen::VectorXd& h,
                                  const std::vector<char>& bounded,bool cholesky=false) {
 using Eigen::MatrixXd;using Eigen::VectorXd;const int q=h.size();
 if(!q) return VectorXd();
 if(q==1) return VectorXd::Constant(1,G(0,0)>0 ? (bounded[0] ? std::max(0.0,h[0]/G(0,0)) : h[0]/G(0,0)) : 0);
 VectorXd d=G.diagonal().cwiseMax(1e-300).cwiseSqrt();
 MatrixXd A=d.cwiseInverse().asDiagonal()*G*d.cwiseInverse().asDiagonal();
 VectorXd b=h.cwiseQuotient(d),x=VectorXd::Zero(q);
 std::vector<char> passive(q,0);
 for(int j=0;j<q;++j) passive[j]=!bounded[j];
 auto solve=[&]() {
   std::vector<int> free;for(int j=0;j<q;++j) if(passive[j]) free.push_back(j);
   MatrixXd H(free.size(),free.size());VectorXd c(free.size()),z=VectorXd::Zero(q);
   for(int a=0;a<(int)free.size();++a) {c[a]=b[free[a]];
     for(int k=0;k<(int)free.size();++k) H(a,k)=A(free[a],free[k]);}
   if(free.size()) {
     VectorXd v;bool solved=false;
     if(cholesky) {
       Eigen::LLT<MatrixXd> factor(H);
       if(factor.info()==Eigen::Success && factor.rcond()>1e-12) {
         v=factor.solve(c);
         // Correct the residual using the same factor, retaining extra precision
         // in its accumulation. No routine QR factorization is needed.
         for(int refinement=0;refinement<2;++refinement) {
           VectorXd residual(c.size());
           for(int i=0;i<c.size();++i) {
             long double value=c[i];
             for(int j=0;j<c.size();++j)value-=(long double)H(i,j)*v[j];
             residual[i]=(double)value;
           }
           VectorXd correction=factor.solve(residual);v+=correction;
           if(correction.lpNorm<Eigen::Infinity>()<=std::numeric_limits<double>::epsilon()*v.lpNorm<Eigen::Infinity>())break;
         }
         const double scale=H.cwiseAbs().rowwise().sum().maxCoeff()*v.lpNorm<Eigen::Infinity>()+
           c.lpNorm<Eigen::Infinity>();
         solved=v.allFinite() && (H*v-c).lpNorm<Eigen::Infinity>()<=1e-11*std::max(1e-300,scale);
       }
     }
     if(!solved)v=H.colPivHouseholderQr().solve(c);
     for(int a=0;a<(int)free.size();++a) z[free[a]]=v[a];}
   return z;
 };
 x=solve();
 const double tol=1e-12*std::max(std::numeric_limits<double>::min(),b.lpNorm<Eigen::Infinity>());
 for(int outer=0;outer<10*q*q;++outer) {
   VectorXd gradient=b-A*x;int enter=-1;double best=tol;
   for(int j=0;j<q;++j) if(!passive[j]&&gradient[j]>best) {best=gradient[j];enter=j;}
   if(enter<0) break;
   passive[enter]=1;
   for(int inner=0;inner<10*q;++inner) {
     VectorXd z=solve();bool positive=true;
     for(int j=0;j<q;++j) if(passive[j]&&bounded[j]&&z[j]<=0) positive=false;
     if(positive) {x=z;break;}
     double fraction=1;
     for(int j=0;j<q;++j) if(passive[j]&&bounded[j]&&z[j]<=0)
       fraction=std::min(fraction,x[j]/std::max(1e-300,x[j]-z[j]));
     x+=fraction*(z-x);
     for(int j=0;j<q;++j) if(passive[j]&&bounded[j]&&x[j]<=1e-14*std::max(std::numeric_limits<double>::min(),x.lpNorm<Eigen::Infinity>()))
       {x[j]=0;passive[j]=0;}
   }
 }
 return x.cwiseQuotient(d);
}
}
#endif
