#ifndef L0GLM_PREDICTOR_GROUPS_H
#define L0GLM_PREDICTOR_GROUPS_H

// Predictor-group passes share the scalar families and storage-aware solvers.
// Only response-channel solves run in OpenMP; R objects stay on the main thread.
struct GroupCacheEntry {MatrixXd gram;VectorXd xy;std::vector<int> active;};
struct GroupCache {std::vector<GroupCacheEntry> channels;};

inline bool group_rowaug(const SparseDesign&,const State&,const std::vector<int>&) {return true;}
inline bool group_rowaug(const DenseDesign& X,const State& s,const std::vector<int>& active) {
 const auto rows=positive_rows(s.W);const double size=(double)rows.size()*active.size();
 if(size<1e6)return false;
 double nonzero=0;
 for(int j:active)for(int i:rows)nonzero+=X(i,j)*std::sqrt(s.W[i])!=0;
 return nonzero<.1*size;
}

inline double group_loglik(const Family& f,const State& s,int kind,const VectorXd& eta) {
 if(kind==7)return NA_REAL;
 long double value=0;
 for(int i=0;i<s.y.size();++i)if(s.w[i]>0) {
  const double y=s.y[i],mu=f.mean(eta[i]),w=s.w[i],phi=s.phi;
  double term=0;
  switch(kind) {
   case 0:term=-.5*(std::log(2*M_PI*phi)-std::log(w)+w*(y-mu)*(y-mu)/phi);break;
   case 1:term=w*R::dpois(std::nearbyint(y),mu,true);break;
   case 2:term=R::dbinom(std::nearbyint(w*y),std::nearbyint(w),mu,true);break;
   case 3:term=w*R::dgamma(y,1/phi,mu*phi,true);break;
   case 4:term=-.5*w*(std::log(2*M_PI*phi*y*y*y)+(y-mu)*(y-mu)/(phi*y*mu*mu));break;
   case 5:term=w*R::dnbinom_mu(y,f.theta,mu,true);break;
   case 6:term=w*((y==0?0:y*std::log(mu))-mu-R::lgammafn(y+1))/phi;break;
  }
  value+=term;
 }
 return (double)value;
}
inline bool group_valid_mu(const Family& f,const VectorXd& mu) {
 for(double v:mu)if(!std::isfinite(v)||(f.kind==2?!(v>0&&v<1):f.kind!=0&&!(v>0)))return false;
 return true;
}
template<class Design> VectorXd group_predict(const Design& X,const VectorXd& b) {
 return predict(X,b,VectorXd::Zero(X.rows()));
}
template<class Design> MatrixXd group_predictions(const Design& X,const MatrixXd& B,
 const std::vector<State>& states,int nthreads) {
 MatrixXd eta(X.rows(),B.cols());
 #ifdef _OPENMP
 #pragma omp parallel for schedule(static) num_threads(nthreads) if(nthreads>1)
 #endif
 for(int k=0;k<B.cols();++k)eta.col(k)=predict(X,VectorXd(B.col(k)),states[k].offset);
 return eta;
}
template<class Design> bool group_cached_ridge(const Design&,State&,const std::vector<int>&,
 const SolverControl&,GroupCacheEntry&,VectorXd&) {return false;}
inline bool group_cached_ridge(const DenseDesign& X,State& s,const std::vector<int>& active,
 const SolverControl& control,GroupCacheEntry& cache,VectorXd& answer) {
 const int p=active.size();std::vector<int> index(X.cols(),-1);
 for(int j=0;j<(int)cache.active.size();++j)index[cache.active[j]]=j;
 bool rebuild=cache.active.empty();for(int a:active)if(index[a]<0)rebuild=true;
 if(rebuild) {
  auto rows=positive_rows(s.w);MatrixXd A=weighted_design(X,rows,active,s.w);
  VectorXd response(rows.size());
  for(int i=0;i<(int)rows.size();++i)response[i]=std::sqrt(s.w[rows[i]])*(s.y[rows[i]]-s.offset[rows[i]]);
  cache.gram=sym_gram(A);cache.xy=right_product(A,MatrixXd(response));cache.active=active;
  for(int j=0;j<p;++j)index[active[j]]=j;
 }
 MatrixXd G(p,p),rhs(p,1),guess(p,1),solution;
 for(int j=0;j<p;++j) {
  const int a=active[j];rhs(j,0)=cache.xy[index[a]]/s.phi+s.lambdas[a]*s.prior[a];guess(j,0)=s.beta[a];
  for(int k=0;k<p;++k)G(j,k)=cache.gram(index[a],index[active[k]])/s.phi;
  G(j,j)+=s.lambdas[a];
 }
 bool success=eigen_coefficients(G,rhs,guess,control,s,solution);
 if(!success&&(s.lambdas.array()==0).count()==p) {
  SolverControl qr;qr.method=3;qr.pivot=1;qr.ordering=1;
  success=eigen_coefficients(G,rhs,guess,qr,s,solution);
 }
 answer=s.beta;
 if(success)for(int j=0;j<p;++j)answer[active[j]]=solution(j,0);
 else s.solver_failed=true;
 return true;
}

template<class Design> List predictor_group_pass(const Design& X,List context,const Family& f,
 const SolverControl& control,int nthreads,SEXP cache_pointer) {
 MatrixXd B=as<MatrixXd>(context["B"]),Y=as<MatrixXd>(context["y"]),
  weights=as<MatrixXd>(context["weights"]),offset=as<MatrixXd>(context["offset"]),
  prior=as<MatrixXd>(context["prior"]),initial_mu=as<MatrixXd>(context["initial_mu"]),
  initial_eta=as<MatrixXd>(context["initial_eta"]);
 LogicalMatrix allowed=clone(as<LogicalMatrix>(context["allowed"]));
 VectorXd lower=as<VectorXd>(context["lower"]),phi=as<VectorXd>(context["phi"]);
 IntegerVector ids=context["ids"],stable=clone(as<IntegerVector>(context["stable"])),included=context["included"];
 LogicalVector estimated=context["estimated"];NumericVector multipliers=context["multiplier"];
 const int n=X.rows(),p=X.cols(),m=B.cols(),G=multipliers.size(),
  maxit=as<int>(context["maxit"]),ll_kind=as<int>(context["loglik_kind"]);
 const double lambda=as<double>(context["lambda"]),thresh=as<double>(context["thresh"]),
  tol=as<double>(context["tol"]),theta_tol=as<double>(context["theta_tol"]);
 const bool joint=as<bool>(context["grouped"]),relaxed=as<bool>(context["relaxed"]),
   auto_rowaug=as<bool>(context["auto_rowaug"]);
 if(n<1||m<1||B.rows()!=p||Y.rows()!=n||Y.cols()!=m||weights.rows()!=n||
    !(weights.cols()==1||weights.cols()==m)||offset.rows()!=n||!(offset.cols()==1||offset.cols()==m)||
    prior.rows()!=p||prior.cols()!=m||lower.size()!=p||allowed.nrow()!=p||allowed.ncol()!=m||
    initial_mu.rows()!=n||initial_mu.cols()!=m||initial_eta.rows()!=n||initial_eta.cols()!=m||
    ids.size()!=p||phi.size()!=m||estimated.size()!=m||stable.size()!=m)
   stop("Incompatible native predictor-group dimensions.");
 for(int j:ids)if(j<1||j>G)stop("Invalid native predictor-group label.");
 for(int j:included)if(j<1||j>p)stop("Invalid native compulsory term.");
 std::vector<char> compulsory(p,0);for(int j:included)compulsory[j-1]=1;
 const int blocks=joint?G:G*m;
 auto block=[&](int j,int k){return ids[j]-1+(joint?0:k*G);};
 GroupCache* cache;
 if(Rf_isNull(cache_pointer)) {
  XPtr<GroupCache> ptr(new GroupCache(),true);cache_pointer=ptr;
  // Protect until the returned list takes ownership of the external pointer.
 }
 Shield<SEXP> protected_cache(cache_pointer);XPtr<GroupCache> ptr(cache_pointer);cache=ptr.get();
 if((int)cache->channels.size()!=m)cache->channels.resize(m);
 std::vector<State> states(m);std::vector<std::string> errors(m);
 for(int k=0;k<m;++k) {
  State& s=states[k];s.y=Y.col(k);s.w=weights.col(weights.cols()==1?0:k);
  s.offset=offset.col(offset.cols()==1?0:k);s.prior=prior.col(k);s.lower=lower;
  s.upper=VectorXd::Constant(p,DBL_MAX);s.scale=VectorXd::Ones(p);
  s.npositive=(s.w.array()>0).count();s.phi=phi[k];
 }
 MatrixXd eta=group_predictions(X,B,states,nthreads);
 std::vector<double> h_ll,h_change,h_noise;std::vector<int> h_selected,h_halvings;
 bool converged=false;int iteration=0,failed_iteration=0,failed_solves=0,inexact=0,fallbacks=0;
 auto loglik=[&](const MatrixXd& e) {
  long double value=0;for(int k=0;k<m;++k)value+=group_loglik(f,states[k],ll_kind,e.col(k));return (double)value;
 };
 auto loss=[&](const MatrixXd& e) {
  long double value=0;for(int k=0;k<m;++k)value+=f.deviance(states[k].y,e.col(k),states[k].w)/states[k].phi;
  return (double)value;
 };
 for(iteration=1;iteration<=maxit;++iteration) {
  if(iteration%8==0)checkUserInterrupt();
  VectorXd norm2=VectorXd::Zero(blocks);
  for(int k=0;k<m;++k)for(int j=0;j<p;++j)if(allowed(j,k)&&!compulsory[j])norm2[block(j,k)]+=B(j,k)*B(j,k);
  MatrixXd D=MatrixXd::Constant(p,m,R_PosInf),proposal=MatrixXd::Zero(p,m);
  for(int k=0;k<m;++k)for(int j=0;j<p;++j)if(compulsory[j])D(j,k)=0;else if(allowed(j,k))
    D(j,k)=lambda==0?0:lambda*multipliers[ids[j]-1]/(4*norm2[block(j,k)]);
  auto channel=[&](int k) {
   try {
    State& s=states[k];s.beta=B.col(k);s.lambdas=D.col(k);s.eta=eta.col(k);
    s.mu.resize(n);for(int i=0;i<n;++i)s.mu[i]=f.mean(s.eta[i]);
    bool valid=group_valid_mu(f,s.mu);for(double e:s.eta)valid=valid&&f.eta_valid(e);
    if(!valid){s.mu=initial_mu.col(k);s.eta=initial_eta.col(k);}
    s.W=VectorXd::Zero(n);s.z=VectorXd::Zero(n);
    for(int i=0;i<n;++i) {
     const double derivative=f.derivative(s.eta[i]);double W,z;
     if(f.link==0){W=s.w[i]/f.variance(s.mu[i]);z=s.y[i]-s.offset[i];}
     else if(f.kind==4&&f.link==1){W=s.w[i]*std::exp(-s.eta[i]);z=s.eta[i]-s.offset[i]-1+s.y[i]/s.mu[i];}
     else {W=s.w[i]*derivative*derivative/f.variance(s.mu[i]);z=s.eta[i]-s.offset[i]+(s.y[i]-s.mu[i])/derivative;}
     if(s.w[i]==0)W=0;
     if(W>0&&std::isfinite(W)&&std::isfinite(z)&&derivative!=0){s.W[i]=W/s.phi;s.z[i]=z;}
    }
    std::vector<int> active;for(int j=0;j<p;++j)if(std::isfinite(s.lambdas[j]))active.push_back(j);
    VectorXd value=VectorXd::Zero(p);
    if(!active.empty()) {
     SolverControl channel_control=control;
     if(auto_rowaug)channel_control.rowaug=group_rowaug(X,s,active);
     const bool cached=f.kind==0&&f.link==0&&!control.external()&&control.method==0&&control.pivot==0&&
       (auto_rowaug||!control.rowaug)&&
       (int)active.size()<=s.npositive&&active.size()<=2000&&active.size()*active.size()*m<=8000000;
     bool solved=cached&&group_cached_ridge(X,s,active,control,cache->channels[k],value);
     if(!solved) {
      bool previous=s.solver_failed;s.solver_failed=false;value=ridge_selected(X,s,active,channel_control,iteration);
      if(s.solver_failed&&control.method==0&&
         std::all_of(active.begin(),active.end(),[&](int j){return s.lambdas[j]==0;})) {
       s.solver_failed=false;SolverControl qr;qr.method=3;qr.pivot=1;qr.ordering=1;
       value=ridge_selected(X,s,active,qr,iteration);
      }
      s.solver_failed=s.solver_failed||previous;
     }
     for(int j=0;j<p;++j)if(!std::isfinite(s.lambdas[j]))value[j]=0;
    }
    for(int j=0;j<p;++j)value[j]=std::max(s.lower[j],value[j]);
    proposal.col(k)=value;
   }catch(const std::exception& e){errors[k]=e.what();}
  };
  if(control.external())for(int k=0;k<m;++k)channel(k);
  else {
   #ifdef _OPENMP
   #pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) if(nthreads>1)
   #endif
   for(int k=0;k<m;++k)channel(k);
  }
  for(const auto& error:errors)if(!error.empty())stop(error);
  auto penalty=[&](const MatrixXd& b) {
   long double value=0;
   for(int k=0;k<m;++k)for(int j=0;j<p;++j)if(std::isfinite(D(j,k)))
    value+=D(j,k)*(b(j,k)-prior(j,k))*(b(j,k)-prior(j,k));
   return (double)value;
  };
  bool changed=false;for(int k=0;k<m;++k)for(int j=0;j<p;++j)changed=changed||((proposal(j,k)!=0)!=(B(j,k)!=0));
  const double baseline=changed?R_PosInf:loss(eta)+penalty(B),
    threshold=std::isfinite(baseline)?baseline+1e-10*std::max(1.,std::abs(baseline)):R_PosInf;
  MatrixXd candidate=proposal,next_eta(n,m);bool accepted=false;int halving=0;
  for(halving=0;halving<=30;++halving) {
   next_eta=group_predictions(X,candidate,states,nthreads);
   const double objective=loss(next_eta)+penalty(candidate);
   if(std::isfinite(objective)&&objective<=threshold){accepted=true;break;}
   candidate=B+std::ldexp(1.,-halving-1)*(proposal-B);
  }
  if(!accepted){failed_iteration=iteration;break;}
  if(!relaxed&&lambda>0) {
   norm2.setZero();for(int k=0;k<m;++k)for(int j=0;j<p;++j)if(allowed(j,k)&&!compulsory[j])
    norm2[block(j,k)]+=candidate(j,k)*candidate(j,k);
   const double largest=norm2.size()?std::sqrt(norm2.maxCoeff()):0;
   MatrixXd cut=candidate;bool cut_changed=false;
   for(int k=0;k<m;++k)for(int j=0;j<p;++j)if(allowed(j,k)&&!compulsory[j]) {
    const double magnitude=std::sqrt(norm2[block(j,k)]);
    if(magnitude<largest*thresh||magnitude<1e-20){cut_changed=cut_changed||cut(j,k)!=0;cut(j,k)=0;}
   }
   MatrixXd cut_eta=next_eta;
   if(cut_changed)cut_eta=group_predictions(X,cut,states,nthreads);
   if(!cut_changed||std::isfinite(loglik(cut_eta))) {
    candidate=cut;next_eta=cut_eta;
    for(int k=0;k<m;++k)for(int j=0;j<p;++j)if(!compulsory[j]) {
     const double magnitude=std::sqrt(norm2[block(j,k)]);
     if(magnitude<largest*thresh||magnitude<1e-20)allowed(j,k)=false;
    }
   }
  }
  double nuisance=0;
  for(int k=0;k<m;++k)if(estimated[k]) {
   bool support_changed=false;int selected=0;
   for(int j=0;j<p;++j){support_changed=support_changed||((candidate(j,k)!=0)!=(B(j,k)!=0));selected+=allowed(j,k);}
   if(support_changed)stable[k]=0;
   if(stable[k]>=3&&!relaxed)continue;
   State& s=states[k];
   if((iteration>=3&&selected<.7*s.npositive)||relaxed) {
    long double rss=0;for(int i=0;i<n;++i)if(s.w[i]>0) {
     const double mu=f.mean(next_eta(i,k)),residual=s.y[i]-mu;rss+=s.w[i]*residual*residual/f.variance(mu);
    }
    double target=s.npositive>selected?(double)rss/(s.npositive-selected):NA_REAL;
    if(std::isfinite(target)&&target>0) {
     stable[k]=!support_changed&&std::abs(target/s.phi-1)<.01?stable[k]+1:0;
     const double updated=std::max(s.phi/2,std::min(2*s.phi,target));
     nuisance=std::max(nuisance,std::abs(std::log(updated/s.phi)));s.phi=updated;
    }
   }
  }
  long double slope_num=0,slope_den=0,pred_num=0,pred_den=0;
  int selected=0;
  for(int j=0;j<p;++j)if(!compulsory[j]) {
   bool present=false;for(int k=0;k<m;++k) {
    const double d=candidate(j,k)-B(j,k);slope_num+=d*d;slope_den+=candidate(j,k)*candidate(j,k);
    present=present||candidate(j,k)!=0;
   }selected+=present;
  }
  std::vector<double> pred_changes(m),pred_sizes(m);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(nthreads) if(nthreads>1)
  #endif
  for(int k=0;k<m;++k) {
   VectorXd difference=group_predict(X,VectorXd(candidate.col(k)-B.col(k))),prediction=group_predict(X,VectorXd(candidate.col(k)));
   pred_changes[k]=difference.squaredNorm();pred_sizes[k]=prediction.squaredNorm();
  }
  for(int k=0;k<m;++k){pred_num+=pred_changes[k];pred_den+=pred_sizes[k];}
  const double change=std::max(std::sqrt((double)slope_num)/std::max(1.,std::sqrt((double)slope_den)),
    std::sqrt((double)pred_num)/std::max(1.,std::sqrt((double)pred_den)));
  B=candidate;eta=next_eta;
  h_ll.push_back(loglik(eta));h_change.push_back(change);h_noise.push_back(nuisance);
  h_selected.push_back(selected);h_halvings.push_back(halving);
  if(change<tol&&nuisance<theta_tol){converged=true;break;}
 }
 for(int k=0;k<m;++k){phi[k]=states[k].phi;failed_solves+=states[k].solver_failed;inexact+=states[k].solver_inexact;fallbacks+=states[k].solver_fallbacks;}
 IntegerVector iters(h_ll.size());for(int i=0;i<iters.size();++i)iters[i]=i+1;
 return List::create(_["B"]=B,_["allowed"]=allowed,_["state"]=List::create(_["phi"]=phi,
   _["estimated"]=estimated,_["stable"]=stable),_["iterations"]=std::min(iteration,maxit),_["converged"]=converged,
   _["history"]=DataFrame::create(_["iteration"]=iters,_["loglik"]=h_ll,_["selected"]=h_selected,
     _["step_halvings"]=h_halvings,_["coefficient_change"]=h_change,_["nuisance_change"]=h_noise),
   _["cache"]=cache_pointer,_["failed_iteration"]=failed_iteration,_["failed_solves"]=failed_solves,
   _["inexact_solves"]=inexact,_["solver_fallbacks"]=fallbacks);
}
#endif
