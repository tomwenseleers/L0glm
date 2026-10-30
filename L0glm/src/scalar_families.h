#ifndef L0GLM_SCALAR_FAMILIES_H
#define L0GLM_SCALAR_FAMILIES_H
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
namespace L0glmScalar {
constexpr double tiny = 1e-15;
struct Terms { double loss, gradient, information; };
struct GammaCache {double phi=-1,a=0,lgamma_a=0,log_a=0;};
inline const GammaCache& gamma_cache(double phi) {
  static thread_local GammaCache c;
  if(phi!=c.phi){c.phi=phi;c.a=1/phi;c.lgamma_a=R::lgammafn(c.a);c.log_a=std::log(c.a);}
  return c;
}
struct CensoredGamma {double logp,q;};
inline CensoredGamma gamma_censor(double mu,double phi,double censor) {
  if(mu==0)return {0,0};
  const auto& c=gamma_cache(phi);const double a=c.a,z=a*censor/mu;
  // Chernoff's upper-tail bound makes this truncation smaller than rounding error.
  if(z>a && -z+a+a*std::log(z/a)<-36) return {0,0};
  if(z/(a+1)<.5) {
    double term=1,sum=1;
    for(int k=1;k<200;++k){term*=z/(a+k);sum+=term;if(term<sum*1e-15)break;}
    return {a*std::log(z)-z-c.lgamma_a-std::log(a)+std::log(sum),a/sum};
  }
  const double lp=R::pgamma(z,a,1.,true,true);
  return {lp,std::exp(a*std::log(z)-z-c.lgamma_a-lp)};
}
inline Terms gamma_scale_terms(double y,double mu,double scale,double censor,bool derivatives=true) {
  if(mu==0)return y==0 ? Terms{0,0,1/tiny} : Terms{R_PosInf,-R_PosInf,R_PosInf};
  const double a=mu/scale;
  if(y>0) return {-R::dgamma(y,a,scale,true),
    derivatives ? (R::digamma(a)+std::log(scale)-std::log(y))/scale : 0,1/(scale*mu)};
  const double z=censor/scale;
  double lp,derivative,information=1/(scale*mu);
  if(z/(a+1)<.5) {
    double term=1,S=1,Sa=0,Saa=0,harmonic=0,harmonic2=0;
    for(int k=1;k<200;++k){term*=z/(a+k);S+=term;
      if(derivatives){const double inv=1/(a+k);harmonic+=inv;harmonic2+=inv*inv;
        Sa-=term*harmonic;Saa+=term*(harmonic*harmonic+harmonic2);}
      if(term<S*1e-15)break;}
    lp=std::min(0.,a*std::log(z)-z-R::lgammafn(a+1)+std::log(S));
    derivative=derivatives ? std::log(z)-R::digamma(a+1)+Sa/S : 0;
    if(derivatives)information=std::max(1e-16,(R::trigamma(a+1)-Saa/S+(Sa/S)*(Sa/S))/(scale*scale));
  } else {
    lp=R::pgamma(z,a,1.,true,true);
    const double step=std::max(1e-7,a*1e-5),lo=std::max(a*.5,a-step),hi=a+step;
    derivative=derivatives ? (R::pgamma(z,hi,1.,true,true)-R::pgamma(z,lo,1.,true,true))/(hi-lo) : 0;
  }
  return {-lp,-derivative/scale,information};
}

// Half minus-two-log-likelihood. Quasi-Poisson uses the scaled Poisson kernel.
inline Terms terms(double y, double mu, int family, double phi, double theta, double censor) {
  if (!std::isfinite(mu) || mu < 0) return {R_PosInf, NA_REAL, NA_REAL};
  if (family == 0) return {.5 * (std::log(2 * M_PI * phi) + (y-mu)*(y-mu)/phi),
    (mu-y)/phi, 1/phi};
  if(family==6)return gamma_scale_terms(y,mu,phi,censor);
  if(family==5) {
    if(phi==1)return terms(y,mu,1,1,theta,censor);
    const double alpha=phi-1;
    if(y==0)return {mu/alpha*std::log(phi),std::log(phi)/alpha,1/(phi*std::max(mu,tiny))};
    if(mu==0)return {R_PosInf,-R_PosInf,R_PosInf};
    const double a=mu/alpha;
    return {-R::dnbinom_mu(y,a,mu,true),
      -(R::digamma(y+a)-R::digamma(a)-std::log(phi))/alpha,1/(phi*mu)};
  }
  if (family == 4 && y == 0) {
    if (mu == 0) return {0, 0, 0};
    const double a=1/phi;const auto c=gamma_censor(mu,phi,censor);
    return {-c.logp,c.q/mu,c.q==0?0:a/(mu*mu)};
  }
  if (mu == 0) {
    if (y == 0) return {0, family == 2 ? 1/phi : 1, 1/tiny};
    return {R_PosInf, -R_PosInf, R_PosInf};
  }
  if (family == 1 || family == 2) {
    const double scale = family == 2 ? phi : 1;
    return {(mu-(y ? y*std::log(mu) : 0)+R::lgammafn(y+1))/scale,
      (1-y/mu)/scale, 1/(scale*mu)};
  }
  if (family == 3) return {
    -R::lgammafn(y+theta)+R::lgammafn(theta)+R::lgammafn(y+1)
      -theta*std::log(theta)-(y ? y*std::log(mu) : 0)+(y+theta)*std::log(mu+theta),
    theta*(mu-y)/(mu*(mu+theta)), theta/(mu*(mu+theta))};
  const auto& cache=gamma_cache(phi);const double a=cache.a;
  return {cache.lgamma_a-a*cache.log_a+a*std::log(mu)
    -(a-1)*std::log(y)+a*y/mu, a*(mu-y)/(mu*mu), a/(mu*mu)};
}

inline Terms derivatives(double y,double mu,int family,double phi,double theta,double censor) {
  if(family==3)return {0,theta*(mu-y)/(mu*(mu+theta)),theta/(mu*(mu+theta))};
  if(family==4 && y>0){const double a=1/phi;return {0,a*(mu-y)/(mu*mu),a/(mu*mu)};}
  if(family==5 && phi>1) {
    const double alpha=phi-1,shape=mu/alpha;
    if(y==0)return {0,std::log1p(alpha)/alpha,0};
    if(y<=64 && y==std::floor(y)) {
      double first=0,second=0;
      for(int k=0;k<(int)y;++k){const double inverse=1/(shape+k);first+=inverse;second+=inverse*inverse;}
      return {0,(std::log1p(alpha)-first)/alpha,second/(alpha*alpha)};
    }
    return {0,-(R::digamma(y+shape)-R::digamma(shape)-std::log(phi))/alpha,
      (R::trigamma(shape)-R::trigamma(shape+y))/(alpha*alpha)};
  }
  if(family==6 && y>0) {
    const double shape=mu/phi;
    return {0,(R::digamma(shape)+std::log(phi)-std::log(y))/phi,
      R::trigamma(shape)/(phi*phi)};
  }
  return terms(y,mu,family,phi,theta,censor);
}

inline double change(double y, double old, double next, int f, double phi, double theta, double censor) {
  if (next == old) return 0;
  if (next < 0 || !std::isfinite(next)) return R_PosInf;
  if (f == 0) return ((y-next)*(y-next)-(y-old)*(y-old))/phi;
  if(f==5 && phi==1)return change(y,old,next,1,1,theta,censor);
  if(f==5 && old>0 && next>0) {
    const double a=old/(phi-1),b=next/(phi-1);
    if(y==0)return 2*(b-a)*std::log(phi);
    if(y<=64 && y==std::floor(y)) {
      double ratio=0;
      for(int k=0;k<(int)y;++k)ratio+=std::log1p((b-a)/(a+k));
      return 2*((b-a)*std::log1p(phi-1)-ratio);
    }
    return 2*(R::lgammafn(b)-R::lgammafn(a)-R::lgammafn(y+b)+R::lgammafn(y+a)+(b-a)*std::log(phi));
  }
  if(f==6 && y>0 && old>0 && next>0) {
    const double a=old/phi,b=next/phi;
    return 2*(R::lgammafn(b)-R::lgammafn(a)+(b-a)*(std::log(phi)-std::log(y)));
  }
  if(f==6)return 2*(gamma_scale_terms(y,next,phi,censor,false).loss-gamma_scale_terms(y,old,phi,censor,false).loss);
  if(f==5)return 2*(terms(y,next,f,phi,theta,censor).loss-terms(y,old,f,phi,theta,censor).loss);
  if (f == 4 && y == 0) return 2*(gamma_censor(old,phi,censor).logp-gamma_censor(next,phi,censor).logp);
  if (next == 0 || old == 0) return 2*(terms(y,next,f,phi,theta,censor).loss-terms(y,old,f,phi,theta,censor).loss);
  const double d = next-old, lr = std::abs(d/old)<.5 ? std::log1p(d/old) : std::log(next)-std::log(old);
  if (f == 1 || f == 2) return 2*(d-y*lr)/(f == 2 ? phi : 1);
  if (f == 3) return 2*((y+theta)*std::log1p(d/(theta+old))-y*lr);
  return 2/phi*(lr+y*(1/next-1/old));
}

}
#endif
