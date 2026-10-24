test_that("native nuisance quadratics and Poisson solves match optimization references", {
qp <- L0glm:::nuisance_quadratic_cpp
set.seed(44)
for(q in 1:7) for(rep in 1:4) {
 A<-matrix(rnorm(30*q),30,q);H<-crossprod(A)+diag(.1,q);h<-rnorm(q)
 bounded<-sample(c(TRUE,FALSE),q,replace=TRUE);lower<-ifelse(bounded,0,-Inf)
 a<-as.numeric(qp(H,h,bounded))
 ref<-optim(rep(0,q),function(b) sum(b*(H%*%b))/2-sum(h*b),
   function(b) as.numeric(H%*%b-h),method="L-BFGS-B",lower=lower,
   control=list(factr=1e3,pgtol=1e-10,maxit=1000))
 expect_equal(sum(a*(H%*%a))/2-sum(h*a),ref$value,tolerance=1e-8)
 g<-as.numeric(H%*%a-h);expect_lt(max(abs(ifelse(bounded&a==0,pmin(g,0),g))),1e-8)
}
tiny<-as.numeric(qp(diag(2),c(1e-17,2e-17),c(TRUE,TRUE)))
expect_lt(max(abs(tiny/c(1e-17,2e-17)-1)),1e-12)
B<-splines::bs(seq(0,1,length.out=120),df=6,intercept=TRUE)
S<-matrix(runif(360),120,3);alpha<-matrix(runif(18,.1,1),6,3)
W<-matrix(runif(360,.2,2),120,3);Y<-matrix(rpois(360,S+B%*%alpha),120,3)
storage.mode(Y)<-"double"
a<-L0glm:::nuisance_poisson_cpp(B,Y,S,W,alpha,60L,1e-8,1L)
b<-L0glm:::nuisance_poisson_cpp(B,Y,S,W,alpha,60L,1e-8,4L)
expect_identical(a,b);expect_lt(max(a$kkt),1e-6)
for(j in 1:3) {
 ref<-optim(alpha[,j],function(v) -sum(W[,j]*dpois(Y[,j],S[,j]+B%*%v,log=TRUE)),
   method="L-BFGS-B",lower=0,control=list(factr=1e3,pgtol=1e-9,maxit=1000))
 expect_equal(-sum(W[,j]*dpois(Y[,j],S[,j]+a$background[,j],log=TRUE)),ref$value,tolerance=1e-7)
}
})
