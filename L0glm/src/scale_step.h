#ifndef L0GLM_SCALE_STEP_H
#define L0GLM_SCALE_STEP_H

// All ordinary estimated-scale families solve the same Pearson equation.
// Guarded selection updates retain the established coefficient-search path.
template<class Mean>
inline double pearson_scale_step(const Family& f,State& s,int k,int iteration,
                                bool support_changed,bool relaxed,Mean mean) {
 if(!s.estimate||s.fixed_scale)return 0;
 if(support_changed){
  if(s.frozen&&!relaxed)++s.restarts;
  s.stable=0;s.frozen=false;
 }
 if(s.npositive<=k||(!relaxed&&(s.frozen||iteration<3||k>=.7*s.npositive)))return 0;
 long double value=0;
 for(int i=0;i<s.y.size();++i)if(s.w[i]>0) {
  const double mu=mean(i),residual=s.y[i]-mu;
  value+=s.w[i]*residual*residual/f.variance(mu);
 }
 const double target=(double)value/(s.npositive-k),old=s.phi;
 if(!(target>0)||!std::isfinite(target))return 0;
 const bool stable=!support_changed&&std::abs(target/old-1)<.01;
 s.stable=stable?s.stable+1:0;
 s.phi=relaxed?target:std::min(2*old,std::max(old/2,target));
 s.frozen=s.stable>=3;++s.updates;
 return std::abs(std::log(s.phi/old));
}

inline double pearson_scale_step(const Family& f,State& s,int k,int iteration,
                                bool support_changed,bool relaxed) {
 return pearson_scale_step(f,s,k,iteration,support_changed,relaxed,[&](int i){return s.mu[i];});
}
#endif
