library(Matrix)

psf = function(x, u, w, h=1) h*exp(((x-u)^2)/(-2*(w^2)))
make_banded_matrix = function(psf, n) { 
                                psf = psf/sum(psf)
                                require(Matrix)
                                bandSparse(n, 
                                k = seq(-(length(psf)-1),0,1), 
                                diag = t(replicate(n, rev(psf))), symm=FALSE)
                           }

x = 1:200
psf_x = psf(0:40, u=20, w=5, h=1)
psf_x[psf_x<1E-5] = 0
plot(psf_x, type="l")
A_psf_x = make_banded_matrix(psf = psf_x, n=length(x))
dim(A_psf_x) # 200 x 200
image(A_psf_x)

y = 1:200
psf_y = psf(0:40, u=20, w=2, h=1)
psf_y[psf_y<1E-5] = 0
plot(psf_y, type="l")
A_psf_y = make_banded_matrix(psf = psf_y, n=length(y))
dim(A_psf_y) # 200 x 200
image(A_psf_y)

image(outer(psf_x, psf_y))

# 2D shifts of 2D PSF across y & x directions calculated as a Kroneker product
# PS using torch there is also https://rdrr.io/cran/GPUmatrix/man/kroneker.html
kronprod = Matrix::kronecker(A_psf_y, A_psf_x) 
dim(kronprod) # 40000 x 40000
100*sum(kronprod!=0) / prod(dim(kronprod)) # 1.6% nonzero
colSums(kronprod)

Y_true = Matrix(0, ncol=length(x), nrow=length(y), sparse=T)
image(Y_true)
nfeatures = 20
xpos = sample.int(n=max(x), size=nfeatures, replace=F)
ypos = sample.int(n=max(y), size=nfeatures, replace=F)
Y_true[cbind(ypos, xpos)] = 10^sapply(1:nfeatures, function (feat) rnorm(1, mean=5, sd=1))
image(Y_true)

Y = as.vector(Y_true) %*% 
