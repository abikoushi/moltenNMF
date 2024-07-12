library(moltenNMF)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS)

center0 <- function(X,B, mu=0){
  i <- attr(X,"assign")
  do.call("rbind", lapply(1:3, function(j)sweep(B[i==j,], 2, colMeans(B[i==j,]) - mu)))  
}


set_data <- function(L, nrow, ncol, mu=0){
  df <- expand.grid(row=1:nrow, col=1:ncol, dep=1:10)
  X <- sparse_onehot(~row+col+dep, data=df)
  B <- matrix(rnorm((nrow+ncol+10)*L,0,1),ncol=L)
  B <- center0(X,B,mu)
  y <- rpois(nrow(X), rowSums(exp(X%*%B)))
  list(X=X, y=y, df=df, trueparam =B)
}
set.seed(599)
dat <- set_data(3, 500, 200, mu = -1)
hist(dat$y, breaks = "scott")
X1 <- with(dat, X[y>0,])
y1 <- with(dat, y[y>0])

attr(X1,"indices") <- attr(dat$X,"indices")

m <- nrow(dat$X)-nrow(X1)
p <- with(dat, colMeans(X[y==0,]))

m/nrow(dat$X)

t1 <- system.time({
  out1 <- moltenNMF::mNMF_vb(dat$y, dat$X,
                              L = 3,
                              iter = 200)
})

t2 <- system.time({
  out2 <- moltenNMF:::mNMF_vb_sp(y1, X1,
                                 L = 3,
                                 N0 = m,
                                 probx = p,
                                 iter = 200)
  out2 <- moltenNMF::mNMF_vb(dat$y, dat$X, 
                             V = out2$shape/out2$rate,
                             L = 3,
                             iter = 1)
})


plot(out1$ELBO, type = "l")
plot(out2$ELBO, type = "l")


f <- moltenNMF::product_m(dat$X, out2$shape/out2$rate)
plot(f, dat$y, pch=".")
abline(0, 1, col="royalblue")

f <- moltenNMF::product_m(dat$X, out1$shape/out1$rate)
plot(f, dat$y, pch=".")
abline(0, 1, col="royalblue")

ind <- apply(cor(dat$trueparam,digamma(out2$shape)-log(out2$rate)),1,which.max)
Bhat <- (digamma(out2$shape)-log(out2$rate))[,ind]
Bhat <- center0(dat$X,Bhat,-1)
plot(dat$trueparam,Bhat,pch=".")
abline(0,1,col="royalblue")

ind <- apply(cor(dat$trueparam,digamma(out1$shape)-log(out1$rate)),1,which.max)
Bhat <- (digamma(out1$shape)-log(out1$rate))[,ind]
Bhat <- center0(dat$X,Bhat,-1)
plot(dat$trueparam,Bhat,pch=".")
abline(0,1,col="royalblue")
