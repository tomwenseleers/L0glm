#ifndef L0GLM_GROUP_H
#define L0GLM_GROUP_H
template<class Design> bool valid_boundary_mean(const Family& f,double mu) {
 if(!std::isfinite(mu))return false;
 if(f.scalar>=0)return (f.scalar==5||f.censor>0)?mu>=0:mu>0;
 if(f.kind==0)return true;
 if(f.kind==2)return mu>0&&mu<1;
 return mu>0;
}
template<class Design> void final_boundary(const Design& X,State& s,const Family& f) {
 VectorXd original=s.beta,previous=s.eta;int removed=0;
 if(s.identity_boundary) {
  bool positive=f.kind==1&&(s.offset.array()>=0).all()&&(s.beta.array()>=0).all()&&L0glm::boundary_nonnegative(X);
  removed=L0glm::boundary_identity(X,s.y,s.w,s.beta,s.eta,s.lower,f.kind==0?0:1,positive);
 }else for(int j=0;j<X.cols();++j)if(s.lower[j]==0&&s.beta[j]>0) {
  double slope=0,old=0,next=0;bool feasible=true,nonempty=false;
  L0glm::boundary_column(X,j,[&](int i,double x) {
   nonempty=true;double eta=s.eta[i]-x*s.beta[j],mu=f.mean(eta);
   if(!f.eta_valid(eta)||!valid_boundary_mean<Design>(f,mu)){feasible=false;return;}
   if(s.w[i]<=0)return;
   double score=(mu-s.y[i])*f.derivative(eta)/f.variance(mu);
   slope+=s.w[i]*x*score;
   old+=s.w[i]*f.loss(s.y[i],f.mean(s.eta[i]),s.eta[i]);
   next+=s.w[i]*f.loss(s.y[i],mu,eta);
  });
  if(!nonempty){s.beta[j]=0;++removed;}
  else if(feasible&&std::isfinite(slope)&&slope>=0&&std::isfinite(next)&&std::isfinite(old)&&next<=old) {
   double b=s.beta[j];s.beta[j]=0;++removed;
   L0glm::boundary_column(X,j,[&](int i,double x){s.eta[i]-=x*b;});
  }
 }
 if(removed) {
  s.eta=predict(X,s.beta,s.offset);
  double value=f.deviance(s.y,s.eta,s.w),old=f.deviance(s.y,previous,s.w);
  if(!std::isfinite(value)||value>old){s.beta=original;s.eta=previous;removed=0;}
 }
 if(removed) {
  for(int i=0;i<X.rows();++i) {
   s.mu[i]=f.mean(s.eta[i]);
   if(f.scalar<0&&f.kind!=0&&f.kind!=2)s.mu[i]=std::max(eps,s.mu[i]);
  }
  initialize_working(s,f,true);
 }
}
template<class Design> void final_state(const Design& X,State& s,const Family& f) {
 if(!s.hard)final_boundary(X,s,f);
 bool any=false;int k=s.offset_df;
 for(int j=0;j<X.cols();++j){any=any||s.beta[j]!=0;k+=s.beta[j]!=0||s.included[j];}
 if(s.estimate&&(!s.hard||s.lambda==0||!any)&&s.npositive>k) {
  long double value=0;
  for(int i=0;i<X.rows();++i)if(s.w[i]>0)value+=s.w[i]*std::pow(s.y[i]-s.mu[i],2)/f.variance(s.mu[i]);
  double phi=(double)value/(s.npositive-k);
  if(std::isfinite(phi)&&phi>0)s.phi=phi;
 }
}
inline VectorXd original_coefficients(const State& s) {
 VectorXd b=s.beta;
 if(s.center_intercept>=0) {
  long double shift=0;for(int j=0;j<b.size();++j)shift+=s.center_shift[j]*b[j];
  b[s.center_intercept]-=(double)shift;
 }
 return b;
}
inline void center_prior(State& s) {
 if(s.center_intercept<0)return;
 s.prior=s.original_prior;long double shift=0;
 for(int j=0;j<s.beta.size();++j)if(s.allowed[j]||s.included[j])shift+=s.center_shift[j]*s.original_prior[j];
 s.prior[s.center_intercept]+=(double)shift;
}
struct PassRecord {int pass,channel,iterations,updates,restarts;bool relaxed,converged;};
inline List pass_result(const State& s) {
 return List::create(_["beta"]=s.beta,_["eta"]=s.eta,_["mu"]=s.mu,_["weights"]=s.Wret,
   _["z"]=s.zret,_["lambdas"]=s.lambdas,_["dispersion"]=s.phi,
   _["iterations"]=s.iterations,_["converged"]=s.converged,_["updates"]=s.updates,_["restarts"]=s.restarts,
   _["failed"]=s.failed,_["solver_failed"]=s.solver_failed,
   _["solver_inexact"]=s.solver_inexact,_["solver_fallbacks"]=s.solver_fallbacks);
}
inline void read_pass(State& s,List result) {
 s.beta=as<VectorXd>(result["beta"]);s.eta=as<VectorXd>(result["eta"]);s.mu=as<VectorXd>(result["mu"]);
 s.Wret=as<VectorXd>(result["weights"]);s.zret=as<VectorXd>(result["z"]);
 s.lambdas=as<VectorXd>(result["lambdas"]);s.phi=as<double>(result["dispersion"]);
 s.iterations=as<int>(result["iterations"]);s.converged=as<bool>(result["converged"]);
 s.updates=as<int>(result["updates"]);s.restarts=as<int>(result["restarts"]);
 s.failed=as<bool>(result["failed"]);
 s.solver_failed=as<bool>(result["solver_failed"]);s.solver_inexact=as<int>(result["solver_inexact"]);
 s.solver_fallbacks=as<int>(result["solver_fallbacks"]);
}
inline List pass_context(const State& s,List base) {
 List context=clone(base);LogicalVector allowed(s.allowed.size());
 for(int j=0;j<allowed.size();++j)allowed[j]=s.allowed[j];
 context["allowed"]=allowed;context["beta"]=s.beta;context["prior.mean"]=s.prior;
 context["auto_beta"]=s.auto_beta;context["auto_phi"]=s.auto_phi;context["dispersion"]=s.phi;
 context["lambda"]=s.lambda;context["hard_threshold"]=s.hard;context["scale_frozen"]=!s.estimate;
 return context;
}
template<class Design> List channel_pass(const Design& X,List contexts,const Family& family,
  const std::vector<double>& probes,const SolverControl& control) {
 auto states=load_states(X,contexts);List result(states.size());
 for(int k=0;k<(int)states.size();++k) {
   prepare_state(X,states[k],family,probes);fit(X,states[k],family,control);
   final_state(X,states[k],family);result[k]=pass_result(states[k]);
 }
 return result;
}
inline List pass_frame(const std::vector<PassRecord>& records) {
 const int n=records.size();CharacterVector stage(n);IntegerVector pass(n),channel(n),iter(n),up(n),rest(n);LogicalVector conv(n);
 for(int j=0;j<n;++j){const auto& z=records[j];stage[j]=z.relaxed?"relaxed":"selection";pass[j]=z.pass;
  channel[j]=z.channel;iter[j]=z.iterations;up[j]=z.updates;rest[j]=z.restarts;conv[j]=z.converged;}
 return DataFrame::create(_["stage"]=stage,_["pass"]=pass,_["channel"]=channel,
  _["iterations"]=iter,_["converged"]=conv,_["dispersion_updates"]=up,_["dispersion_restarts"]=rest);
}
template<class Design> List model(const Design& X,List contexts,const Family& family,
 const std::vector<double>& probes,int nthreads,int max_passes,bool relaxed,bool nonnegative,
 const SolverControl& control,SEXP executor=R_NilValue) {
 const int p=X.cols(),m=contexts.size(),n=X.rows();
 std::vector<State> states=load_states(X,contexts);
 std::vector<std::vector<char>> fixed_allowed(m);std::vector<VectorXd> initial_eta(m),initial_mu(m);
 for(int j=0;j<m;++j){fixed_allowed[j]=states[j].allowed;initial_eta[j]=states[j].eta;initial_mu[j]=states[j].mu;center_prior(states[j]);}
 std::vector<char> active(p,0),last_active(p,0);
 for(int j=0;j<p;++j)for(const auto& s:states)active[j]=active[j]||s.allowed[j]||s.included[j];
 std::vector<PassRecord> history;std::vector<std::string> errors(m);
 int failed_steps=0,failed_solves=0,inexact_solves=0,solver_fallbacks=0;
 MatrixXd original=MatrixXd::Zero(p,m);
 auto pool=[&]() {
  for(int k=0;k<m;++k)original.col(k)=original_coefficients(states[k]);
  std::vector<double> seeds(p,0);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(nthreads) if(p>10000)
  #endif
  for(int j=0;j<p;++j) {
   double average=0,magnitude=0,maximum=-1,direction=0;bool present=false;
   for(int k=0;k<m;++k) {
    double b=original(j,k);average+=b;magnitude+=std::abs(b);present=present||b!=0;
    if(std::abs(b)>maximum){maximum=std::abs(b);direction=b>0?1:b<0?-1:0;}
   }
   average/=m;magnitude/=m;
   active[j]=present||states[0].included[j];
   seeds[j]=nonnegative?(average==0&&present?magnitude:average):
     (average>0?1:average<0?-1:direction)*magnitude;
  }
  for(int k=0;k<m;++k) {
   State& s=states[k];s.beta.setZero();
   for(int j=0;j<p;++j){s.allowed[j]=active[j]&&fixed_allowed[k][j];if(s.allowed[j]||s.included[j])s.beta[j]=seeds[j];}
   if(s.center_intercept>=0) {
    long double shift=0;for(int j=0;j<p;++j)shift+=s.center_shift[j]*s.beta[j];
    s.beta[s.center_intercept]+=(double)shift;
   }
   center_prior(s);s.auto_beta=false;s.auto_phi=false;
  }
 };
 auto run=[&](int pass,bool refit) {
  last_active=active;
   for(int k=0;k<m;++k) {
   State& s=states[k];
   s.eta=initial_eta[k];s.mu=initial_mu[k];s.hard=!refit;s.lambda=refit?1e-30:s.selection_lambda;
   s.stable=s.iterations=s.updates=s.restarts=0;
    s.converged=s.failed=s.solver_failed=false;s.frozen=!s.estimate;s.old_prediction.resize(0);
    s.solver_inexact=s.solver_fallbacks=0;
   }
   if(executor!=R_NilValue) {
     List jobs(m);for(int k=0;k<m;++k)jobs[k]=pass_context(states[k],contexts[k]);
     Function execute(executor);List results=execute(jobs);
     for(int k=0;k<m;++k)read_pass(states[k],results[k]);
   }else {
     #ifdef _OPENMP
     #pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) if(nthreads>1)
     #endif
     for(int k=0;k<m;++k)try {prepare_state(X,states[k],family,probes);}
       catch(const std::exception& e){errors[k]=e.what();}
     for(const auto& e:errors)if(!e.empty())stop(e);
     prefill_shared_gram(X,states,family);
     auto solve_channel=[&](int k) {
       try {fit(X,states[k],family,control);final_state(X,states[k],family);}
       catch(const std::exception& e){errors[k]=e.what();}
     };
     if(control.external())for(int k=0;k<m;++k)solve_channel(k);
     else {
       #ifdef _OPENMP
       #pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) if(nthreads>1)
       #endif
       for(int k=0;k<m;++k)solve_channel(k);
     }
     for(const auto& e:errors)if(!e.empty())stop(e);
   }
  for(int k=0;k<m;++k) {
   const auto& s=states[k];history.push_back({pass,k+1,s.iterations,s.updates,s.restarts,refit,s.converged});
    failed_steps+=s.failed;failed_solves+=s.solver_failed;
    inexact_solves+=s.solver_inexact;solver_fallbacks+=s.solver_fallbacks;
  }
 };
 for(int pass=1;pass<=max_passes;++pass) {
  if(pass>1)pool();
  if(pass>1&&!std::any_of(active.begin(),active.end(),[](char x){return x;}))break;
  run(pass,false);
 }
 if(relaxed) {
  pool();if(std::any_of(active.begin(),active.end(),[](char x){return x;}))run(1,true);
 }
 MatrixXd coefficients(p,m),penalties=MatrixXd::Zero(p,m),weights(n,m),z(n,m);NumericVector dispersion(m);
 for(int k=0;k<m;++k) {
  coefficients.col(k)=original_coefficients(states[k]);weights.col(k)=states[k].Wret;z.col(k)=states[k].zret;
  dispersion[k]=states[k].phi;
  for(int j=0;j<p;++j)if(last_active[j])penalties(j,k)=states[k].lambdas[j];
 }
 return List::create(_["coefficients"]=coefficients,_["lambdas"]=penalties,_["weights"]=weights,
  _["z"]=z,_["dispersion"]=dispersion,_["convergence"]=pass_frame(history),
  _["native_stats"]=List::create(_["failed_steps"]=failed_steps,_["failed_solves"]=failed_solves,
     _["inexact_solves"]=inexact_solves,_["solver_fallbacks"]=solver_fallbacks));
}
#endif
