  static MatrixXd contrasts(int C) {
    MatrixXd Q=MatrixXd::Zero(C,C-1);
    for(int j=0;j<C-1;++j){double norm=std::sqrt((j+1.)*(j+2.));for(int k=0;k<=j;++k)Q(k,j)=-1/norm;Q(j+1,j)=(j+1)/norm;}
    return Q;
  }
  double multinomial_loglik(const MatrixXd& B) {
    const auto& eta=prediction(B);long double value=0;
    for(int i=0;i<n;++i)if(weights(i,0)>0) {
      double maximum=eta.row(i).maxCoeff();if(!std::isfinite(maximum))return R_NegInf;
      double total=Y.row(i).sum(),normalizer=(eta.row(i).array()-maximum).exp().sum();
      double term=R::lgammafn(total+1)-total*std::log(normalizer);
      for(int c=0;c<m;++c)term+=Y(i,c)*(eta(i,c)-maximum)-R::lgammafn(Y(i,c)+1);
      value+=weights(i,0)*term;
    }
    return (double)value;
  }
  struct ContrastBlock {MatrixXd T;VectorXd start,prior;int offset=0;bool free=false;};
  MatrixXd multinomial_working(const MatrixXd& B,const MatrixXd& D) {
    std::vector<int> columns,rows;
    for(int j=0;j<p;++j)if((D.row(j).array().isFinite()).count()>1)columns.push_back(j);
    MatrixXd answer=MatrixXd::Zero(p,m);if(columns.empty())return answer;
    for(int i=0;i<n;++i)if(weights(i,0)*Y.row(i).sum()>0)rows.push_back(i);
    const int N=rows.size(),P=columns.size(),K=m-1;MatrixXd Q=contrasts(m);
    Matrix A0=weighted_design(X,rows,columns,VectorXd::Ones(n));
    const auto& eta=prediction(B);MatrixXd probability(N,m),score(N,m),linear(N,K),d(P,m);
    VectorXd weight(N);bool shared=true;
    for(int i=0;i<N;++i) {
      int row=rows[i];double maximum=eta.row(row).maxCoeff(),total=Y.row(row).sum();
      probability.row(i)=(eta.row(row).array()-maximum).exp();probability.row(i)/=probability.row(i).sum();
      score.row(i)=weights(row,0)*(Y.row(row)-total*probability.row(i));weight[i]=weights(row,0)*total;
      Eigen::RowVectorXd predictor=eta.row(row);
      // A one-column offset is shared across all classes and cancels in Q.
      if(offset.cols()==1)predictor.array()-=offset(row,0);else predictor-=offset.row(row);
      linear.row(i)=predictor*Q;
    }
    for(int j=0;j<P;++j)for(int c=0;c<m;++c) {
      d(j,c)=D(columns[j],c);shared=shared&&std::isfinite(d(j,c))&&d(j,c)==D(columns[j],0);
    }
    std::vector<ContrastBlock> blocks(P);VectorXd start,center,penalty;int dimension=0;
    SparseMatrix<double> map;
    if(shared) {
      MatrixXd b(P,m),a(P,m);for(int j=0;j<P;++j){b.row(j)=B.row(columns[j]);a.row(j)=prior.row(columns[j]);}
      MatrixXd seed=b*Q,location=a*Q;start=Eigen::Map<VectorXd>(seed.data(),seed.size());center=Eigen::Map<VectorXd>(location.data(),location.size());
      dimension=P*K;penalty.resize(dimension);for(int c=0;c<K;++c)penalty.segment(c*P,P)=d.col(0);
    }else {
      std::vector<Eigen::Triplet<double>> entries;std::vector<double> seeds,centers,penalties;
      for(int j=0;j<P;++j) {
        std::vector<int> keep;for(int c=0;c<m;++c)if(std::isfinite(d(j,c)))keep.push_back(c);
        const int size=keep.size(),rank=size-1;auto& block=blocks[j];block.offset=dimension;block.free=true;
        VectorXd diagonal(size),b(size),a(size);for(int c=0;c<size;++c){diagonal[c]=d(j,keep[c]);b[c]=B(columns[j],keep[c]);a[c]=prior(columns[j],keep[c]);block.free=block.free&&diagonal[c]==0;}
        MatrixXd q=contrasts(size),T;VectorXd seed;
        MatrixXd information=q.transpose()*diagonal.asDiagonal()*q;
        Eigen::LLT<MatrixXd> factor(information);
        if(block.free||factor.info()==Eigen::Success) {
          MatrixXd R=block.free?MatrixXd::Identity(rank,rank):MatrixXd(factor.matrixU());
          T=q*R.triangularView<Eigen::Upper>().solve(MatrixXd::Identity(rank,rank));seed=R*(q.transpose()*b);
        }else {
          VectorXd root=diagonal.cwiseSqrt(),v=root.cwiseInverse();v/=v.maxCoeff();v/=v.norm();v[0]+=1;
          MatrixXd U=(MatrixXd::Identity(size,size)-2*v*v.transpose()/v.squaredNorm()).rightCols(rank);
          T=root.cwiseInverse().asDiagonal()*U;seed=U.transpose()*root.cwiseProduct(b);
        }
        block.T=MatrixXd::Zero(m,rank);for(int c=0;c<size;++c)block.T.row(keep[c])=T.row(c);
        VectorXd location=block.free?VectorXd::Zero(rank):VectorXd(T.transpose()*diagonal.cwiseProduct(a));
        for(int c=0;c<rank;++c){seeds.push_back(seed[c]);centers.push_back(location[c]);penalties.push_back(block.free?0:1);}
        for(int c=0;c<m;++c)for(int k=0;k<rank;++k)if(block.T(c,k)!=0)entries.emplace_back(j+c*P,dimension+k,block.T(c,k));
        dimension+=rank;
      }
      map.resize(P*m,dimension);map.setFromTriplets(entries.begin(),entries.end());
      start=Eigen::Map<VectorXd>(seeds.data(),seeds.size());center=Eigen::Map<VectorXd>(centers.data(),centers.size());penalty=Eigen::Map<VectorXd>(penalties.data(),penalties.size());
    }
    if(std::min(dimension,N*K)>max_working)stop("Multinomial working system exceeds control.family$max.working.");
    VectorXd coefficient;State diagnostics;
    if(dimension<=N*K) {
      L0glm::require_workspace(48.*dimension*dimension,control.workspace,"Multinomial Hessian");
      MatrixXd H=MatrixXd::Zero(dimension,dimension);VectorXd gradient;
      auto cross=[&](const VectorXd& w) {
        Matrix weighted=A0;
        if constexpr(std::is_same_v<Matrix,SparseMatrix<double>>) {
          for(int j=0;j<weighted.outerSize();++j)for(typename Matrix::InnerIterator it(weighted,j);it;++it)it.valueRef()*=w[it.row()];
        }else for(int j=0;j<P;++j)weighted.col(j).array()*=w.array();
        return MatrixXd(A0.transpose()*weighted);
      };
      if(shared) {
        MatrixXd mean=probability*Q;std::vector<std::pair<int,int>> pairs;
        for(int j=0;j<K;++j)for(int k=0;k<=j;++k)pairs.emplace_back(j,k);
        std::vector<std::string> errors(pairs.size());
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,1) num_threads(threads) if(threads>1)
        #endif
        for(int index=0;index<(int)pairs.size();++index) {
          try {
          int j=pairs[index].first,k=pairs[index].second;
          VectorXd value=probability*(Q.col(j).cwiseProduct(Q.col(k)));
          value.array()-=mean.col(j).array()*mean.col(k).array();value.array()*=weight.array();
          MatrixXd block=cross(value);H.block(j*P,k*P,P,P)=block;H.block(k*P,j*P,P,P)=block.transpose();
          }catch(const std::exception& error){errors[index]=error.what();}
        }
        for(const auto& error:errors)if(!error.empty())stop(error);
        MatrixXd g=A0.transpose()*(score*Q);gradient=Eigen::Map<VectorXd>(g.data(),g.size());
      }else {
        for(int j=0;j<m;++j)for(int k=0;k<m;++k) {
          VectorXd value=(weight.array()*probability.col(j).array()*((j==k?1.:0.)-probability.col(k).array())).matrix();
          MatrixXd block=cross(value);SparseMatrix<double> left=map.middleRows(j*P,P),right=map.middleRows(k*P,P);
          H.noalias()+=left.transpose()*(block*right);
        }
        MatrixXd g=A0.transpose()*score;gradient=map.transpose()*Eigen::Map<VectorXd>(g.data(),g.size());
      }
      MatrixXd rhs=H*start+gradient+penalty.cwiseProduct(center),solution;H.diagonal()+=penalty;
      if(!eigen_coefficients(H,rhs,MatrixXd(start),control,diagnostics,solution))stop("Native multinomial quadratic solver failed.");
      coefficient=solution.col(0);
    }else {
      // Build the exact C-1 whitened observation operator, retaining sparse X.
      if constexpr(std::is_same_v<Matrix,MatrixXd>)L0glm::require_workspace(48.*N*K*dimension,control.workspace,"Multinomial working design");
      std::vector<MatrixXd> roots(N);VectorXd response(N*K);
      for(int i=0;i<N;++i) {
        VectorXd v=probability.row(i).transpose().cwiseMax(1e-12);v/=v.sum();VectorXd mean=Q.transpose()*v;
        MatrixXd V=weight[i]*(Q.transpose()*v.asDiagonal()*Q-mean*mean.transpose());Eigen::LLT<MatrixXd> factor(V);
        if(factor.info()!=Eigen::Success)stop("Native multinomial contrast information is not positive definite.");
        MatrixXd R=factor.matrixU();roots[i]=R;
        VectorXd z=R*linear.row(i).transpose()+R.transpose().triangularView<Eigen::Lower>().solve(VectorXd(Q.transpose()*score.row(i).transpose()));
        for(int c=0;c<K;++c)response[i+c*N]=z[c];
      }
      Matrix design;std::vector<Eigen::Triplet<double>> entries;
      if constexpr(std::is_same_v<Matrix,MatrixXd>)design=MatrixXd::Zero(N*K,dimension);
      auto add=[&](int i,int j,double value) {
        MatrixXd T=shared?roots[i]:MatrixXd(roots[i]*Q.transpose()*blocks[j].T);
        for(int a=0;a<K;++a)for(int b=0;b<T.cols();++b) {
          int column=shared?j+b*P:blocks[j].offset+b;double number=value*T(a,b);
          if constexpr(std::is_same_v<Matrix,MatrixXd>)design(i+a*N,column)=number;
          else if(number!=0)entries.emplace_back(i+a*N,column,number);
        }
      };
      if constexpr(std::is_same_v<Matrix,SparseMatrix<double>>) {
        for(int j=0;j<P;++j)for(typename Matrix::InnerIterator it(A0,j);it;++it)add(it.row(),j,it.value());
        design.resize(N*K,dimension);design.setFromTriplets(entries.begin(),entries.end());
      }else for(int j=0;j<P;++j)for(int i=0;i<N;++i)add(i,j,A0(i,j));
      if(!selected_dual(design,response,penalty,center,control,diagnostics,coefficient))stop("Native multinomial Woodbury solver failed.");
    }
    if(shared) {
      Eigen::Map<MatrixXd> beta(coefficient.data(),P,K);MatrixXd value=beta*Q.transpose();
      for(int j=0;j<P;++j)answer.row(columns[j])=value.row(j);
    }else for(int j=0;j<P;++j)answer.row(columns[j])=(blocks[j].T*coefficient.segment(blocks[j].offset,blocks[j].T.cols())).transpose();
    return answer;
  }
