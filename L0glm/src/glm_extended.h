#ifndef L0GLM_EXTENDED_NATIVE_H
#define L0GLM_EXTENDED_NATIVE_H

// Preparation and R callbacks stay on the R thread. Numerical channel work uses
// only native data and the same storage-aware ridge kernels as ordinary GLMs.
template<class Design> struct ExtendedModel {
  using Matrix=std::conditional_t<std::is_same_v<Design,DenseDesign>,MatrixXd,SparseMatrix<double>>;
  const Design& X;SolverControl control;Family family;
  int n,p,m,threads,kind,max_working;
  Eigen::Map<MatrixXd> Y,weights,offset,prior;
  VectorXd lower,theta,noise_df;std::vector<int> stable,estimated,ids;
  std::vector<double> multiplier;std::vector<char> compulsory;
  double theta_min,theta_max;bool fixed_noise,predictor_groups;
  std::vector<State> states;std::vector<L0glm::RidgeWorkspace<Matrix>> workspace;
  MatrixXd cached_B,cached_eta;bool prediction_ready=false;
  static Eigen::Map<MatrixXd> view(List context,const char* name) {
    NumericMatrix matrix=context[name];return Eigen::Map<MatrixXd>(matrix.begin(),matrix.nrow(),matrix.ncol());
  }
  ExtendedModel(const Design& x,List c,const Family& f,const SolverControl& solver,int cores):
    X(x),control(solver),family(f),n(x.rows()),p(x.cols()),m(NumericMatrix(c["B"]).ncol()),
    threads(std::max(1,std::min(cores,m))),kind(as<int>(c["kind"])),max_working(as<int>(c["max_working"])),
    Y(view(c,"y")),weights(view(c,"weights")),offset(view(c,"offset")),
    prior(view(c,"prior")),lower(as<VectorXd>(c["lower"])),theta(as<VectorXd>(c["theta"])),noise_df(as<VectorXd>(c["noise_df"])),
    stable(as<std::vector<int>>(c["stable"])),estimated(as<std::vector<int>>(c["estimated"])),
    ids(as<std::vector<int>>(c["ids"])),multiplier(as<std::vector<double>>(c["multiplier"])),
    theta_min(as<double>(c["theta_min"])),theta_max(as<double>(c["theta_max"])),
    fixed_noise(as<bool>(c["fixed_noise"])),predictor_groups(as<bool>(c["predictor_groups"])),
    states(m),workspace(threads) {
    if(Y.rows()!=n||Y.cols()!=m||prior.rows()!=p||prior.cols()!=m||lower.size()!=p||
      weights.rows()!=n||offset.rows()!=n||!(weights.cols()==1||weights.cols()==m)||
      !(offset.cols()==1||offset.cols()==m)||ids.size()!=(size_t)p)stop("Invalid native extended-family dimensions.");
    compulsory.assign(p,0);for(int j:as<std::vector<int>>(c["included"]))compulsory.at(j-1)=1;
    for(int& id:ids)--id;
    for(int k=0;k<m;++k) {
      auto& s=states[k];s.y=Y.col(k);s.w=weights.col(weights.cols()==1?0:k);
      s.offset=offset.col(offset.cols()==1?0:k);s.prior=prior.col(k);s.lower=lower;
      s.upper=VectorXd::Constant(p,DBL_MAX);s.scale=VectorXd::Ones(p);s.phi=1;s.workspace=control.workspace;
    }
    if(kind==2)cox_setup(c,cores);
  }
  const MatrixXd& prediction(const MatrixXd& B) {
    if(prediction_ready&&cached_B.rows()==B.rows()&&(cached_B.array()==B.array()).all())return cached_eta;
    cached_B=B;cached_eta=group_predictions(X,B,states,threads);prediction_ready=true;return cached_eta;
  }
  double linkfun(double mu)const {
    switch(family.link){case 0:return mu;case 1:return std::log(mu);case 2:return std::sqrt(mu);
      case 3:return 1/mu;case 8:return 1/(mu*mu);}
    return NA_REAL;
  }
  #include "glm_multinomial.h"
  #include "glm_cox.h"
  double loglik(const MatrixXd& B) {
    if(kind==1)return multinomial_loglik(B);
    if(kind==2)return cox_stats(B).loglik;
    const MatrixXd& eta=prediction(B);std::vector<double> values(m,0);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads) if(threads>1)
    #endif
    for(int k=0;k<m;++k) {
      auto f=family;f.theta=theta[k];bool valid=true;
      for(int i=0;i<n;++i)if(states[k].w[i]>0) {
        double mu=f.mean(eta(i,k));if(!std::isfinite(mu)||mu<=0){valid=false;break;}
      }
      values[k]=valid?group_loglik(f,states[k],5,eta.col(k)):R_NegInf;
    }
    long double value=0;for(double v:values)value+=v;return (double)value;
  }
  MatrixXd working(const MatrixXd& B,const MatrixXd& D) {
    if(kind==1)return multinomial_working(B,D);
    if(kind==2)return cox_working(B,D);
    const MatrixXd& eta=prediction(B);MatrixXd answer=MatrixXd::Zero(p,m);std::vector<std::string> errors(m);
    auto run=[&](int k,int worker) {
      try {
        State& s=states[k];auto f=family;f.theta=theta[k];s.beta=B.col(k);s.lambdas=D.col(k);s.eta=eta.col(k);
        s.mu.resize(n);bool valid=true;
        for(int i=0;i<n;++i){s.mu[i]=f.mean(s.eta[i]);valid=valid&&std::isfinite(s.mu[i])&&s.mu[i]>0;}
        if(!valid)for(int i=0;i<n;++i){s.mu[i]=s.y[i]+.1;s.eta[i]=linkfun(s.mu[i]);}
        initialize_working(s,f,false);std::vector<int> active;
        for(int j=0;j<p;++j)if(std::isfinite(D(j,k)))active.push_back(j);
        VectorXd value=VectorXd::Zero(p);
        if(!active.empty()) {
          value=ridge_selected(X,s,active,control,1,&workspace[worker]);
          if(s.solver_failed)throw std::runtime_error("Native negative-binomial ridge solve failed.");
          for(int j=0;j<p;++j)if(!std::isfinite(D(j,k)))value[j]=0;
        }
        answer.col(k)=value.cwiseMax(lower);
      }catch(const std::exception& e){errors[k]=e.what();}
    };
    if(control.external())for(int k=0;k<m;++k)run(k,0);
    else {
      #ifdef _OPENMP
      #pragma omp parallel for schedule(static) num_threads(threads) if(threads>1)
      #endif
      for(int k=0;k<m;++k) {
        int worker=0;
        #ifdef _OPENMP
        worker=omp_get_thread_num();
        #endif
        run(k,worker);
      }
    }
    for(const auto& error:errors)if(!error.empty())stop(error);
    return answer;
  }
  MatrixXd projected(const MatrixXd& B,const MatrixXd& D) {
    const MatrixXd& eta=prediction(B);MatrixXd answer=B;
    for(int k=0;k<m;++k) {
      auto f=family;f.theta=theta[k];const auto& s=states[k];VectorXd gradient_rows(n),information(n);
      for(int i=0;i<n;++i) {
        double mu=f.mean(eta(i,k)),derivative=f.derivative(eta(i,k));
        if(!std::isfinite(mu)||mu<=0)return B;
        double w=s.w[i]/f.variance(mu);gradient_rows[i]=w*derivative*(mu-s.y[i]);information[i]=w*derivative*derivative;
      }
      VectorXd gradient=X.transpose()*gradient_rows;
      for(int j=0;j<p;++j)if(std::isfinite(D(j,k))) {
        double diagonal=D(j,k),g=gradient[j]+D(j,k)*(B(j,k)-prior(j,k));
        if constexpr(std::is_same_v<Matrix,SparseMatrix<double>>) {
          for(typename Design::InnerIterator it(X,j);it;++it)diagonal+=information[it.row()]*it.value()*it.value();
        }else for(int i=0;i<n;++i)diagonal+=information[i]*X(i,j)*X(i,j);
        if(diagonal>0&&std::isfinite(diagonal)&&std::isfinite(g))answer(j,k)=std::max(lower[j],B(j,k)-g/diagonal);
      }
    }
    return answer;
  }
  double update_noise(const MatrixXd& B,const MatrixXd& previous,int iteration,bool refit) {
    if(kind!=0||fixed_noise)return 0;
    const MatrixXd& eta=prediction(B);std::vector<double> changes(m,0);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads) if(threads>1)
    #endif
    for(int k=0;k<m;++k)if(estimated[k]) {
      const auto& s=states[k];int selected=0,count=0;double weight=0;bool changed=false,valid=true;
      for(int j=0;j<p;++j){selected+=B(j,k)!=0||compulsory[j];changed=changed||((B(j,k)!=0)!=(previous(j,k)!=0));}
      if(noise_df[k]>=0)selected=noise_df[k];
      if(changed)stable[k]=0;
      VectorXd mu(n);for(int i=0;i<n;++i) {
        mu[i]=family.mean(eta(i,k));if(s.w[i]>0){++count;weight+=s.w[i];valid=valid&&mu[i]>0&&std::isfinite(mu[i]);}
      }
      if(!valid||count<=selected||!(refit||(stable[k]<3&&iteration>=3&&selected<.7*count)))continue;
      const double target=weight/count*(count-selected);
      auto equation=[&](double log_theta) {
        long double total=0;double t=std::exp(log_theta);
        for(int i=0;i<n;++i)if(s.w[i]>0)total+=s.w[i]*std::pow(s.y[i]-mu[i],2)/(mu[i]+mu[i]*mu[i]/t);
        return (double)total/target-1;
      };
      double lo=std::log(theta_min),hi=std::log(theta_max),a=equation(lo),b=equation(hi),target_theta;
      if(!std::isfinite(a)||!std::isfinite(b))continue;
      if(b<=0)target_theta=theta_max;else if(a>=0)target_theta=theta_min;
      else {for(int it=0;it<64&&hi-lo>1e-8;++it){double mid=(lo+hi)/2;if(equation(mid)>0)hi=mid;else lo=mid;}target_theta=std::exp((lo+hi)/2);}
      const double old=theta[k];stable[k]=!changed&&std::abs(target_theta/old-1)<.01?stable[k]+1:0;
      theta[k]=refit?target_theta:std::min(2*old,std::max(old/2,target_theta));changes[k]=std::abs(std::log(theta[k]/old));
    }
    return *std::max_element(changes.begin(),changes.end());
  }
  List state()const {
    LogicalVector estimate(m);for(int k=0;k<m;++k)estimate[k]=estimated[k]!=0;
    return List::create(_["theta"]=theta,_["estimated"]=estimate,_["stable"]=stable);
  }
};

template<class Design> List extended_pass(const Design& X,List context,const Family& family,
    const SolverControl& control,int nthreads) {
  ExtendedModel<Design> model(X,context,family,control,nthreads);
  if(as<bool>(context["initialize"])) {
    MatrixXd initial=model.cox_initial();return List::create(_["B"]=initial,_["diagnostic"]=model.cox_info(true));
  }
  MatrixXd B=as<MatrixXd>(context["B"]);LogicalMatrix r_allowed=context["allowed"];
  Eigen::ArrayXXi allowed(B.rows(),B.cols());for(int k=0;k<B.cols();++k)for(int j=0;j<B.rows();++j)allowed(j,k)=r_allowed(j,k);
  const int p=B.rows(),m=B.cols(),G=model.multiplier.size(),maxit=as<int>(context["maxit"]);
  const double lambda=as<double>(context["lambda"]),thresh=as<double>(context["thresh"]),
    tol=as<double>(context["tol"]),theta_tol=as<double>(context["theta_tol"]);
  const bool grouped=as<bool>(context["grouped"]),refit=as<bool>(context["refit"]),pilot=as<bool>(context["pilot"]);
  std::vector<double> h_ll,h_change,h_noise,h_residual;std::vector<int> h_selected,h_halvings,h_inner;
  std::vector<std::string> h_backend;
  bool converged=false;int iteration=0,failed_iteration=0;
  const int blocks=grouped?G:G*m;
  auto block=[&](int j,int k){return model.ids[j]+(grouped?0:k*G);};
  for(iteration=1;iteration<=maxit;++iteration) {
    if(iteration%8==0)checkUserInterrupt();
    VectorXd norm=VectorXd::Zero(blocks);
    for(int k=0;k<m;++k)for(int j=0;j<p;++j)if(allowed(j,k)&&!model.compulsory[j])norm[block(j,k)]+=B(j,k)*B(j,k);
    MatrixXd D=MatrixXd::Constant(p,m,R_PosInf);
    for(int k=0;k<m;++k)for(int j=0;j<p;++j) {
      if(model.compulsory[j])D(j,k)=0;
      else if(allowed(j,k))D(j,k)=pilot?1e-3:lambda==0?0:lambda*model.multiplier[model.ids[j]]/(4*norm[block(j,k)]);
    }
    if(pilot)D=as<MatrixXd>(context["D"]);
    MatrixXd proposal=model.working(B,D);
    if(pilot) {
      List result=List::create(_["B"]=proposal);
      if(model.kind==2)result["diagnostic"]=model.cox_info();
      return result;
    }
    auto objective=[&](const MatrixXd& b) {
      long double ridge=0;for(int k=0;k<m;++k)for(int j=0;j<p;++j)if(std::isfinite(D(j,k)))
        ridge+=D(j,k)*std::pow(b(j,k)-model.prior(j,k),2);
      return -2*model.loglik(b)+(double)ridge;
    };
    bool support_changed=((proposal.array()!=0)!=(B.array()!=0)).any();
    const bool boundary_refit=refit&&(model.lower.array()==0).any();
    double baseline=boundary_refit||!support_changed?objective(B):R_PosInf;
    MatrixXd candidate;int halvings=0;double candidate_value=R_PosInf;
    auto search=[&](const MatrixXd& trial,int limit,double rtol) {
      for(int h=0;h<=limit;++h) {
        candidate=h==0?trial:B+std::ldexp(1.,-h)*(trial-B);
        try {candidate_value=objective(candidate);}
        catch(const std::exception&){candidate_value=R_PosInf;}
        if(std::isfinite(candidate_value)&&candidate_value<=baseline+rtol*std::max(1.,std::abs(baseline))) {halvings+=h;return true;}
      }
      halvings+=limit+1;return false;
    };
    bool accepted=search(proposal,boundary_refit?10:30,boundary_refit?1e-12:1e-10);
    if(boundary_refit&&std::isfinite(baseline)&&(!accepted||candidate_value>=baseline-1e-12*std::max(1.,std::abs(baseline)))) {
      accepted=search(model.projected(B,D),18,1e-12);if(!accepted){candidate=B;accepted=true;}
    }
    if(!accepted){failed_iteration=iteration;break;}
    if(!refit&&lambda>0) {
      norm.setZero();for(int k=0;k<m;++k)for(int j=0;j<p;++j)if(allowed(j,k)&&!model.compulsory[j])norm[block(j,k)]+=candidate(j,k)*candidate(j,k);
      VectorXd largest=VectorXd::Zero(m);
      for(int k=0;k<m;++k)for(int j=0;j<p;++j)if(!model.compulsory[j])largest[k]=std::max(largest[k],std::sqrt(norm[block(j,k)]));
      if(grouped||model.predictor_groups)largest.setConstant(largest.maxCoeff());
      MatrixXd cut=candidate;auto cut_allowed=allowed;
      for(int k=0;k<m;++k)for(int j=0;j<p;++j)if(allowed(j,k)&&!model.compulsory[j]) {
        double magnitude=std::sqrt(norm[block(j,k)]);
        if(magnitude<largest[k]*thresh||magnitude<1e-20){cut(j,k)=0;cut_allowed(j,k)=0;}
      }
      if(model.kind==1)for(int j=0;j<p;++j) {
        int count=0;double mean=0;for(int k=0;k<m;++k)if(cut_allowed(j,k)){++count;mean+=cut(j,k);}
        if(count<2){cut.row(j).setZero();cut_allowed.row(j).setZero();}
        else for(int k=0;k<m;++k)if(cut_allowed(j,k))cut(j,k)-=mean/count;
      }
      if(std::isfinite(model.loglik(cut))){candidate=cut;allowed=cut_allowed;}
    }
    const double noise=model.update_noise(candidate,B,iteration,refit||lambda==0);
    long double num=0,den=0;int selected=0;
    for(int j=0;j<p;++j)if(!model.compulsory[j]) {
      bool present=false;for(int k=0;k<m;++k){double d=candidate(j,k)-B(j,k);num+=d*d;den+=candidate(j,k)*candidate(j,k);present=present||candidate(j,k)!=0;}selected+=present;
    }
    std::vector<double> pred_change(m),pred_size(m);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(model.threads) if(model.threads>1)
    #endif
    for(int k=0;k<m;++k) {
      pred_change[k]=group_predict(X,VectorXd(candidate.col(k)-B.col(k))).squaredNorm();
      pred_size[k]=group_predict(X,VectorXd(candidate.col(k))).squaredNorm();
    }
    long double pn=0,pd=0;for(int k=0;k<m;++k){pn+=pred_change[k];pd+=pred_size[k];}
    double change=std::max(std::sqrt((double)num)/std::max(1.,std::sqrt((double)den)),
      std::sqrt((double)pn)/std::max(1.,std::sqrt((double)pd)));
    B=candidate;h_ll.push_back(model.loglik(B));h_change.push_back(change);h_noise.push_back(noise);
    h_selected.push_back(selected);h_halvings.push_back(halvings);
    if(model.kind==2){h_backend.push_back(model.cox_diagnostic.backend);h_inner.push_back(model.cox_diagnostic.iterations);h_residual.push_back(model.cox_diagnostic.residual);}
    if(change<tol&&noise<theta_tol){converged=true;break;}
  }
  LogicalMatrix output_allowed(p,m);for(int k=0;k<m;++k)for(int j=0;j<p;++j)output_allowed(j,k)=allowed(j,k);
  IntegerVector iters(h_ll.size());for(int i=0;i<iters.size();++i)iters[i]=i+1;
  DataFrame history=DataFrame::create(_["iteration"]=iters,_["loglik"]=h_ll,_["selected"]=h_selected,
    _["step_halvings"]=h_halvings,_["coefficient_change"]=h_change,_["nuisance_change"]=h_noise);
  if(model.kind==2){history["backend"]=h_backend;history["inner_iterations"]=h_inner;history["inner_residual"]=h_residual;}
  return List::create(_["B"]=B,_["allowed"]=output_allowed,_["state"]=model.state(),
    _["converged"]=converged,_["iterations"]=std::min(iteration,maxit),_["failed_iteration"]=failed_iteration,_["history"]=history);
}
#endif
