library(moltenNMF)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS)
library(readr)

center0 <- function(X,B, mu=0){
  i <- attr(X,"assign")
  do.call("rbind", lapply(1:3, function(j)sweep(B[i==j,], 2, colMeans(B[i==j,]) - mu)))  
}


set_data <- function(L, nrow, ncol, ndep=10, mu=0){
  df <- expand.grid(row=1:nrow, col=1:ncol, dep=1:ndep)
  X <- sparse_onehot(~row+col+dep, data=df)
  B <- matrix(rnorm((nrow+ncol+10)*L,0,1),ncol=L)
  B <- center0(X,B,mu)
  y <- rpois(nrow(X), rowSums(exp(X%*%B)))
  list(X=X, y=y, df=df, trueparam = B)
}

set.seed(99); dat <- set_data(3, 500, 200, mu = -0.5)
m <- nrow(dat$X)-nrow(X1)
p <- with(dat, colMeans(X[y==0,]))

cat(m, "\n", file = "testdat.csv", append = FALSE)
cat(p, "\n", file = "testdat.csv", append = TRUE)
cat(attr(dat$X, "indices"), "\n", file = "testdat.csv",
    append = TRUE)
cbind(dat$df, y=dat$y) %>% 
  dplyr::filter(y > 0) %>% 
  write_csv(file = "testdat.csv", append = TRUE, 
            col_names = TRUE)



ind <- scan("testdat.csv", nmax = 4, what = integer())
dat <- read_csv(con, n_max = 10, skip=1)


X1 <- with(dat, X[y>0,])
y1 <- with(dat, y[y>0])

attr(X1,"indices") <- attr(dat$X,"indices")



m/nrow(dat$X)

t1 <- system.time({
  out1 <- moltenNMF::mNMF_vb(dat$y, dat$X,
                              L = 3,
                              iter = 100)
})
 
t2 <- system.time({
  out2 <- moltenNMF:::mNMF_vb_sp(y1, X1,
                                 L = 3,
                                 N0 = m,
                                 probx = p,
                                 iter = 97)
  elbo2 <-out2$ELBO
  out2 <- moltenNMF::mNMF_vb(dat$y, dat$X, 
                             V = out2$shape/out2$rate,
                             L = 3,
                             iter = 3)
  elbo2 <- c(elbo2,out2$ELBO)
})

t1
t2
plot(out1$ELBO, type = "l")
lines(elbo2, type = "l")


f <- moltenNMF::product_m(dat$X, out2$shape/out2$rate)
plot(f, dat$y, pch=".", main = "skip-zeros")
abline(0, 1, col="royalblue")

f <- moltenNMF::product_m(dat$X, out1$shape/out1$rate)
plot(f, dat$y, pch=".", main = "ordinary")
abline(0, 1, col="royalblue")

Bhat <- digamma(out2$shape)-log(out2$rate)
Bhat <- center0(dat$X,Bhat,-0.5)
ind <- apply(cor(dat$trueparam,Bhat),1,which.max)
Bhat <- Bhat[,ind]
dim(Bhat)

plot(dat$trueparam,Bhat,pch=".")
abline(0,1,col="royalblue")

ind <- apply(cor(dat$trueparam,digamma(out1$shape)-log(out1$rate)),1,which.max)
Bhat <- (digamma(out1$shape)-log(out1$rate))[,ind]
Bhat <- center0(dat$X,Bhat,-0.5)
plot(dat$trueparam,Bhat,pch=".")
abline(0,1,col="royalblue")
