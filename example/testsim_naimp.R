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

Y=dat$Y
ind_na = sample.int(110*100, 1000)
Y[ind_na] <- NA

out <- moltenNMF:::NMF2D_vb(Y, rank = 2, iter = 500)
plot(out$ELBO, type="l")
V = moltenNMF:::meanV_2D(out)
plot(log1p(V[[1]]%*%t(V[[2]])), as.matrix(log1p(dat$Y)), pch=1, cex=0.5, col=rgb(0,0,0,0.2))
abline(0, 1, col="royalblue",lty=2)

fit = V[[1]]%*%t(V[[2]])
Yimp = Y
plot(dat$Y[ind_na], fit[ind_na])
abline(0,1)
image(as.matrix(log1p(Y)), col=hcl.colors(12))
