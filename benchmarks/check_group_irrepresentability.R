# Fixed-design group-LASSO diagnostic for the website GC/MS simulation.
# See Bach (2008), JMLR 9:1179-1225, equations (4)-(5).
library(Matrix)
z <- readRDS("codex_sandbox/website_large_examples/multi-comparison.rds")
X <- z$data$X
B <- as.matrix(z$beta_true)
w <- 1/(rowMeans(z$data$y)+0.1)
A <- Diagonal(x=sqrt(w)) %*% X
A <- A %*% Diagonal(x=1/sqrt(colSums(A^2)))
S <- which(rowSums(B^2)>0)
N <- setdiff(seq_len(ncol(A)),S)
AS <- A[,S,drop=FALSE]
G <- crossprod(AS)
U <- B[S,,drop=FALSE]/sqrt(rowSums(B[S,,drop=FALSE]^2))
V <- solve(G,U)
residual <- max(abs(G %*% V-U))/max(abs(U))
stopifnot(is.finite(residual),residual<1e-7)
C <- as.matrix(crossprod(A[,N,drop=FALSE],AS %*% V))
norms <- sqrt(rowSums(C^2))
# Positive-part scores describe the inactive-row KKT bound. They do not include
# the tangent-cone adjustment for zero entries inside active nonnegative groups.
positive_norms <- sqrt(rowSums(pmax(C,0)^2))
answer <- data.frame(condition=c("standard_group_condition","positive_part_diagnostic"),
  maximum=c(max(norms),max(positive_norms)),
  rows_exceeding_one=c(sum(norms>1+1e-7),sum(positive_norms>1+1e-7)),
  relative_solve_residual=residual,n=nrow(A),p=ncol(A),channels=ncol(B),
  true_rows=length(S))
print(answer)
write.csv(answer,"benchmarks/plots/website/multichannel-irrepresentability.csv",
  row.names=FALSE)
