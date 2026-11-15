  struct CoxValues {VectorXd e,S0,q,score;double loglik=R_NegInf;};
  struct CoxSolve {VectorXd value;int iterations=0;double residual=0;bool success=false;std::string backend="pcg",space="";int dimension=0,halvings=0;};
  VectorXd cox_d,cox_event;std::vector<int> cox_ends,cox_group;
  double cox_tol=1e-10,cox_cutoff=1024;int cox_maxit=2000,cox_cores=1;std::string cox_mode="auto";
  CoxSolve cox_diagnostic;
  CoxValues cox_values;MatrixXd cox_B,cox_means;bool cox_ready=false;
  Matrix cox_weighted;std::vector<int> cox_weighted_active;
  mutable VectorXd cox_u,cox_risk,cox_back,cox_response;
  void cox_setup(List context,int cores) {
    List cache=context["cox_cache"];
    cox_d=as<VectorXd>(cache["d"]);cox_event=as<VectorXd>(cache["event"]);
    cox_ends=as<std::vector<int>>(cache["ends"]);cox_group=as<std::vector<int>>(cache["group"]);
    for(int& i:cox_ends)--i;
    for(int& i:cox_group)--i;
    cox_tol=as<double>(context["cox_tol"]);cox_maxit=as<int>(context["cox_maxit"]);
    cox_cutoff=as<double>(context["cox_cutoff"]);cox_mode=as<std::string>(context["cox_mode"]);cox_cores=std::max(1,cores);
  }
  const CoxValues& cox_stats(const MatrixXd& B) {
    if(cox_ready&&cox_B.rows()==B.rows()&&(cox_B.array()==B.array()).all())return cox_values;
    cox_B=B;cox_ready=true;
    const auto& eta=prediction(B);auto& out=cox_values;const int groups=cox_d.size();out.loglik=R_NegInf;
    out.e.setZero(n);out.S0.setZero(groups);out.q.setZero(n);out.score.setZero(n);
    double shift=R_NegInf;for(int i=0;i<n;++i){if(!std::isfinite(eta(i,0)))return out;if(states[0].w[i]>0)shift=std::max(shift,eta(i,0));}
    double sum=0;int group=0;
    for(int i=0;i<n;++i) {
      if(states[0].w[i]>0)out.e[i]=states[0].w[i]*std::exp(eta(i,0)-shift);
      sum+=out.e[i];if(i==cox_ends[group]){out.S0[group]=sum;++group;}
    }
    VectorXd back(groups);sum=0;
    for(int g=groups-1;g>=0;--g) {
      if(cox_d[g]>0&&out.S0[g]<=0)stop("Cox risk-set sums underflowed; use less extreme starts.");
      if(out.S0[g]>0)sum+=cox_d[g]/out.S0[g];
      back[g]=sum;
    }
    long double ll=0;for(int i=0;i<n;++i) {
      out.q[i]=out.e[i]*back[cox_group[i]];out.score[i]=cox_event[i]-out.q[i];
      if(cox_event[i]>0)ll+=cox_event[i]*(eta(i,0)-shift-std::log(out.S0[cox_group[i]]));
    }
    out.loglik=(double)ll;return out;
  }
  template<class Operator> CoxSolve cox_cg(Operator apply,const VectorXd& rhs,const VectorXd& diagonal)const {
    CoxSolve out;out.dimension=rhs.size();VectorXd scale=diagonal.cwiseSqrt(),b=rhs.cwiseQuotient(scale);
    auto scaled=[&](const VectorXd& v){return VectorXd(apply(VectorXd(v.cwiseQuotient(scale))).cwiseQuotient(scale));};
    VectorXd value=VectorXd::Zero(rhs.size()),r=b,direction=r;double norm=b.norm(),rr=r.squaredNorm();
    if(norm==0){out.value=value;out.success=true;return out;}
    for(int iteration=1;iteration<=cox_maxit;++iteration) {
      out.iterations=iteration;VectorXd product=scaled(direction);double curvature=direction.dot(product);
      if(!(curvature>0)||!std::isfinite(curvature))break;
      double alpha=rr/curvature;value+=alpha*direction;r-=alpha*product;double next=r.squaredNorm();
      if(std::sqrt(next)<=cox_tol*norm) {
        r=b-scaled(value);next=r.squaredNorm();
        if(std::sqrt(next)<=cox_tol*norm){out.success=true;break;}
        direction=r;rr=next;continue;
      }
      direction=r+(next/rr)*direction;rr=next;
    }
    out.value=value.cwiseQuotient(scale);out.residual=(b-scaled(value)).norm()/norm;
    out.success=out.success&&out.value.allFinite();return out;
  }
  template<class Selected> VectorXd cox_product(const Selected& A,const CoxValues& stats,
      const VectorXd& d,const VectorXd& v)const {
    cox_u.noalias()=A*v;cox_risk.resize(cox_d.size());auto& u=cox_u;auto& risk=cox_risk;double sum=0;int g=0;
    for(int i=0;i<n;++i){sum+=stats.e[i]*u[i];if(i==cox_ends[g]){risk[g]=sum;++g;}}
    cox_back.resize(cox_d.size());auto& back=cox_back;sum=0;
    for(int j=cox_d.size()-1;j>=0;--j) {
      if(stats.S0[j]>0)sum+=cox_d[j]*risk[j]/(stats.S0[j]*stats.S0[j]);
      back[j]=sum;
    }
    cox_response.resize(n);auto& response=cox_response;
    for(int i=0;i<n;++i)response[i]=stats.q[i]*u[i]-stats.e[i]*back[cox_group[i]];
    return A.transpose()*response+d.cwiseProduct(v);
  }
  template<class Selected> MatrixXd cox_information(const Selected& A,const CoxValues& stats) {
    const int size=A.cols();std::vector<int> events;for(int j=0;j<cox_d.size();++j)if(cox_d[j]>0)events.push_back(j);
    cox_means.resize(events.size(),size);auto& means=cox_means;auto& weighted=cox_weighted;
    if constexpr(std::is_same_v<Matrix,SparseMatrix<double>>) {
      for(int j=0;j<weighted.outerSize();++j)for(int at=weighted.outerIndexPtr()[j];at<weighted.outerIndexPtr()[j+1];++at)
        weighted.valuePtr()[at]=A.valuePtr()[at]*stats.q[weighted.innerIndexPtr()[at]];
    }else for(int j=0;j<size;++j)weighted.col(j).array()=A.col(j).array()*stats.q.array();
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(cox_cores) if(cox_cores>1&&size>32)
    #endif
    for(int j=0;j<size;++j) {
      double sum=0;int k=0;
      for(int i=0;i<n&&k<(int)events.size();++i) {
        sum+=stats.e[i]*A.coeff(i,j);
        if(i==cox_ends[events[k]]){int g=events[k];means(k,j)=sum/stats.S0[g]*std::sqrt(cox_d[g]);++k;}
      }
    }
    MatrixXd information;
    if constexpr(std::is_same_v<Matrix,MatrixXd>)information=right_product(A,weighted);
    else information=MatrixXd(A.transpose()*weighted);
    return information-sym_gram(means);
  }
  MatrixXd cox_working(const MatrixXd& B,const MatrixXd& D) {
    MatrixXd answer=MatrixXd::Zero(p,1);std::vector<int> active,rows(n);std::iota(rows.begin(),rows.end(),0);
    for(int j=0;j<p;++j)if(std::isfinite(D(j,0)))active.push_back(j);
    int size=active.size();if(!size){cox_diagnostic=CoxSolve();cox_diagnostic.backend="empty";cox_diagnostic.success=true;return answer;}
    const auto& stats=cox_stats(B);VectorXd d(size),beta(size),center(size);
    for(int j=0;j<size;++j){int a=active[j];d[j]=D(a,0);beta[j]=B(a,0);center[j]=prior(a,0);}
    auto& work=workspace[0];work.update(X,rows,active,VectorXd::Ones(n),VectorXd::Ones(size),VectorXd::Zero(n));const auto& A=work.A;
    bool unrestricted=(lower.array()<-1e100).all();int events=(cox_d.array()>0).count();
    if(cox_mode=="pcg"&&!unrestricted)stop("Cox cox='pcg' requires unrestricted coefficients.");
    bool direct_ok=size<=n&&size<=max_working&&(double)size*events<=1e7&&48.*size*size+8.*size*events<=control.workspace;
    bool direct=cox_mode=="direct"||!unrestricted||(cox_mode=="auto"&&size<=cox_cutoff&&direct_ok);
    VectorXd proposal;
    if(!direct) {
      VectorXd rhs=A.transpose()*stats.score-d.cwiseProduct(beta-center),diagonal=d;
      for(int j=0;j<size;++j) {
        if constexpr(std::is_same_v<Matrix,SparseMatrix<double>>) {
          for(typename Matrix::InnerIterator it(A,j);it;++it)diagonal[j]+=stats.q[it.row()]*it.value()*it.value();
        }else diagonal[j]+=A.col(j).array().square().matrix().dot(stats.q);
        if(diagonal[j]<=0)diagonal[j]=1;
      }
      cox_diagnostic=cox_cg([&](const VectorXd& v){return cox_product(A,stats,d,v);},rhs,diagonal);
      if(cox_diagnostic.success)proposal=beta+cox_diagnostic.value;
      else if(cox_mode=="auto"&&direct_ok) {
        warning("Cox PCG did not converge; falling back to the direct solver.");direct=true;
      }else stop("Cox PCG did not converge; increase cox.cg.maxit or use cox='direct'.");
    }
    if(direct) {
      if(std::min(size,n)>max_working)stop("Cox working system exceeds control.family$max.working; screen predictors or use cox='pcg'.");
      if(size>n) {
        L0glm::require_workspace(48.*n*n+48.*n*size,control.workspace,"Cox wide working design");
        VectorXd back(cox_d.size());double sum=0;
        for(int g=cox_d.size()-1;g>=0;--g){if(stats.S0[g]>0)sum+=cox_d[g]/(stats.S0[g]*stats.S0[g]);back[g]=sum;}
        MatrixXd W(n,n);for(int j=0;j<n;++j)for(int i=0;i<n;++i)
          W(i,j)=(i==j?stats.q[i]:0)-stats.e[i]*stats.e[j]*back[std::max(cox_group[i],cox_group[j])];
        Eigen::SelfAdjointEigenSolver<MatrixXd> factor(W);
        if(factor.info()!=Eigen::Success)stop("Cox wide information factorization failed.");
        std::vector<int> keep;double largest=factor.eigenvalues().maxCoeff();
        for(int j=0;j<n;++j)if(factor.eigenvalues()[j]>largest*1e-13)keep.push_back(j);
        MatrixXd root(keep.size(),n);VectorXd response(keep.size());
        const auto& eta=prediction(B);VectorXd linear=eta.col(0)-states[0].offset;
        for(int j=0;j<(int)keep.size();++j) {
          int index=keep[j];double scale=std::sqrt(factor.eigenvalues()[index]);
          root.row(j)=scale*factor.eigenvectors().col(index).transpose();
          response[j]=root.row(j).dot(linear)+factor.eigenvectors().col(index).dot(stats.score)/scale;
        }
        MatrixXd design=root*A;State diagnostic;
        if(!selected_dual(design,response,d,center,control,diagnostic,proposal))stop("Cox wide ridge solve failed.");
        cox_diagnostic=CoxSolve();cox_diagnostic.backend="direct_wide";cox_diagnostic.success=true;
        for(int j=0;j<size;++j)answer(active[j],0)=std::max(lower[active[j]],proposal[j]);
        return answer;
      }
      L0glm::require_workspace(48.*size*size+8.*size*events,control.workspace,"Cox direct information");
      if(cox_weighted_active!=active){cox_weighted=A;cox_weighted_active=active;}
      MatrixXd H=cox_information(A,stats),rhs=H*beta+A.transpose()*stats.score+d.cwiseProduct(center),solution;
      H.diagonal()+=d;State diagnostic;
      if(!eigen_coefficients(H,rhs,MatrixXd(beta),control,diagnostic,solution))stop("Cox quadratic solver failed.");
      proposal=solution.col(0);bool fallback=!cox_diagnostic.success&&cox_diagnostic.iterations>0;
      cox_diagnostic.backend=fallback?"direct_fallback":"llt";cox_diagnostic.success=true;
      if(!fallback){cox_diagnostic.iterations=0;cox_diagnostic.residual=NA_REAL;}
    }
    for(int j=0;j<size;++j)answer(active[j],0)=std::max(lower[active[j]],proposal[j]);
    return answer;
  }
  MatrixXd cox_initial() {
    MatrixXd zero=MatrixXd::Zero(p,1);const auto& stats=cox_stats(zero);std::vector<int> rows,columns(p);std::iota(columns.begin(),columns.end(),0);
    for(int i=0;i<n;++i)if(stats.q[i]>0)rows.push_back(i);
    bool wide=p>(int)rows.size();const double ridge=1e-3;
    VectorXd z=VectorXd::Zero(n);for(int i:rows)z[i]=stats.score[i]/stats.q[i];
    auto& work=workspace[0];work.update(X,rows,columns,stats.q,VectorXd::Ones(p),z);const auto& A=work.A;
    VectorXd beta;
    if(wide) {
      VectorXd diagonal=VectorXd::Ones(rows.size());
      if constexpr(std::is_same_v<Matrix,SparseMatrix<double>>) {
        for(int j=0;j<A.outerSize();++j)for(typename Matrix::InnerIterator it(A,j);it;++it)diagonal[it.row()]+=it.value()*it.value()/ridge;
      }else diagonal.array()+=A.array().square().rowwise().sum()/ridge;
      cox_diagnostic=cox_cg([&](const VectorXd& v){return VectorXd(v+A*(A.transpose()*v)/ridge);},work.response,diagonal);
      beta=A.transpose()*cox_diagnostic.value/ridge;
    }else {
      VectorXd diagonal=VectorXd::Constant(p,ridge);for(int j=0;j<p;++j)diagonal[j]+=A.col(j).squaredNorm();
      VectorXd rhs=A.transpose()*work.response;
      cox_diagnostic=cox_cg([&](const VectorXd& v){return VectorXd(A.transpose()*(A*v)+ridge*v);},rhs,diagonal);
      beta=cox_diagnostic.value;
    }
    if(!cox_diagnostic.success)stop("The matrix-free Cox ridge initializer did not converge; increase cox.cg.maxit.");
    cox_diagnostic.space=wide?"observations":"coefficients";MatrixXd proposal=beta.cwiseMax(lower),candidate;
    const double baseline=-2*stats.loglik;bool accepted=false;
    for(int h=0;h<=30;++h) {
      candidate=std::ldexp(1.,-h)*proposal;double objective=-2*cox_stats(candidate).loglik+ridge*candidate.squaredNorm();
      if(std::isfinite(objective)&&objective<=baseline){accepted=true;cox_diagnostic.halvings=h;break;}
    }
    if(!accepted)stop("No descent step for the Cox ridge initializer.");
    return candidate;
  }
  List cox_info(bool initial=false)const {
    List out=List::create(_["backend"]=cox_diagnostic.backend,_["iterations"]=cox_diagnostic.iterations,
      _["residual"]=cox_diagnostic.residual,_["converged"]=cox_diagnostic.success);
    if(initial){out["space"]=cox_diagnostic.space;out["dimension"]=cox_diagnostic.dimension;out["step_halvings"]=cox_diagnostic.halvings;}
    return out;
  }
