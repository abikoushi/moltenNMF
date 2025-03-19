library(moltenNMF)
library(Matrix)
library(dplyr)
library(ggplot2)

set_data_mf <- function(L, nrow, ncol, mu=0){
  W <- matrix(rnorm(nrow*L,0,1),ncol=L)
  H <- matrix(rnorm(L*ncol,0,1),nrow=L)
  W <- sweep(W,1,rowMeans(W)-mu)
  H <- sweep(H,2,rowMeans(H)-mu)
  Y <- matrix(rpois(nrow*ncol, exp(W%*%H)), nrow, ncol)
  Y <- as(Y, "TsparseMatrix")
  list(Y=Y, trueW=W, trueH=H)
}

dat <- set_data_mf(2,110, 130)
hist(dat$Y@x)
#lr = 1000/nnzero(dat$Y)
out <- moltenNMF:::NMF2D_svb(dat$Y, rank = 2,
                             n_epochs = 50, n_baches = as.integer(2000),
                             lr_param = c(1.5,0.6), lr_type = "power")
head(out$shape[[1]])
head(out$shape[[2]])
out$rate
plot(out$ELBO, type="l")
V = moltenNMF:::meanV_2D(out)
plot(V[[1]]%*%t(V[[2]]), as.matrix(dat$Y), pch=1, cex=0.5, col=rgb(0,0,0,0.2))
abline(0, 1, col="royalblue",lty=2)