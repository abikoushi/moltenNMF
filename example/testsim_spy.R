library(Matrix)
library(dplyr)
library(ggplot2)
library(moltenNMF)

set_data_mf <- function(L, nrow, ncol, mu=0){
  W <- matrix(rnorm(nrow*L,0,1),ncol=L)
  H <- matrix(rnorm(L*ncol,0,1),nrow=L)
  W <- sweep(W,1,rowMeans(W)-mu)
  H <- sweep(H,2,colMeans(H)-mu)
  Y <- matrix(rpois(nrow*ncol, exp(W)%*%exp(H)), nrow, ncol)
  Y <- as(Y, "TsparseMatrix")
  list(Y=Y, trueW=W, trueH=H)
}

dat <- set_data_mf(2, 502, 11)
X <- moltenNMF::sparse_onehot(~row+col, data=expand.grid(row=1:502, col=1:11))

bm = bench::mark({
  out_d <- moltenNMF:::mNMF_vb.default(as.integer(dat$Y), X = X, L = 2, iter=1000)
},iterations = 1)

plot(out_d$ELBO[-1], type = "l")
V <- out_d$shape/out_d$rate
f_d <- moltenNMF::product_m(X, V)
length(f_d)
length(dat$Y)
plot(as.matrix(dat$Y), f_d,  pch=1, col=rgb(0,0,0,0.2), xlab="fitted", ylab="obsereved")
abline(0, 1, col="grey", lty=2)

#with skip-zeros
y = as.integer(dat$Y)
wch = which(y>0)
Y = sparseVector(y[wch], wch, length = length(y))
bm2 = bench::mark({
  out <- moltenNMF:::mNMF_vb.default(Y, X = X, L = 2, iter=1000)
}, iterations = 1)

plot(out_d$ELBO[-1], type = "l")
lines(out$ELBO[-1], type = "l", lty=2)
bm$time
bm2$time
bm$mem_alloc
bm2$mem_alloc

Y1 = y[wch]
X1 = slice_rows(X, wch)
X0 = X[-wch,]
# X1@p
# max(X1@p)+max(X0@p)
# dim(X)
# X0@i <- sample(X1@i, size = length(X0@i))
# 
# all(colSums(X[-wch,])==diff(X0@p))
# 

system.time({
  out2 = moltenNMF:::mNMF_vb_sp(y = Y1, X = X1, xp0 = X0@p, L=2,
                                iter=1000,
                                a=0.5, b=0.01,
                                V=NULL,
                                display_progress=TRUE,
                                indices=NULL)  
})

# plot(out_d$ELBO[-1], type = "l")
# plot(out2$ELBO[-1], type = "l", col="red")

V <- out2$shape/out2$rate
f2 <- moltenNMF::product_m(X, V)
plot(c(as.matrix(dat$Y)), f2,  pch=1, col=rgb(0,0,0,0.1), xlab="fitted", ylab="obsereved",cex=0.5)
points(c(as.matrix(dat$Y)), f_d,  pch=2, col=rgb(0,0.5,1,0.1), xlab="fitted", ylab="obsereved", cex=0.5)
abline(0, 1, col="grey", lty=2)

mean((c(as.matrix(dat$Y))-f2)^2)
mean((c(as.matrix(dat$Y))-f_d)^2)

# probX0 = colMeans(X[-wch,])
# N0 = nrow(X)-length(wch)
# 
# out2 = moltenNMF:::mNMF_vb_sp(y = Y1, X = X1, N = N0, probX0 = probX0, L=2, 
#                        iter=1000,
#                        a=0.5, b=0.01,
#                        V=NULL,
#                        display_progress=TRUE,
#                        indices=NULL)
# plot(out2$ELBO, type = "l")

V <- out2$shape/out2$rate
f2 <- moltenNMF::product_m(X, V)
plot(as.matrix(dat$Y), f2,  pch=1, col=rgb(0,0,0,0.2), xlab="fitted", ylab="obsereved")
points(as.matrix(dat$Y), f_d,  pch=2, col=rgb(0,0,0.5,0.2), xlab="fitted", ylab="obsereved")
abline(0, 1, col="grey", lty=2)


####
##SVB omit 0s
y = as.integer(dat$Y)
wch = which(y>0)
Y1 = y[wch]
X1 = slice_rows(X, wch)

probX0 = colMeans(X[-wch,])
N0 = nrow(X)-length(wch)

system.time({
  out_s <- moltenNMF:::mNMF_svb_sp(Y1, X = X1,
                                   N = nrow(X), L = 2,
                                   n_epochs = 1, 
                                   n_batches = 2000,
                                   lr_param = c(15,0.9), 
                                   lr_type = "exponential",
                                   display_progress = TRUE)
})
# user  system elapsed 
# 11.477   0.057  11.654 
plot(out_s$ELBO, type="l")
V <- out_s$shape/out_s$rate
f <- moltenNMF::product_m(X, V)
plot(as.matrix(dat$Y), f,  pch=2, col=rgb(0,0,0,0.5))
points(as.matrix(dat$Y), f_d,  pch=1, col=rgb(1,0.5,0,0.5))
abline(0, 1, col="grey", lty=2)

# plot(as.numeric(log1p(dat$Y)), log1p(f_d),  pch=1, col=rgb(1,0.5,0,0.1))
# points(as.numeric(log1p(dat$Y)), log1p(f),  pch=2, col=rgb(0,0.5,1,0.1))
# abline(0, 1, col="grey", lty=2)

##########
y = as.integer(dat$Y)
wch = which(y>0)
Y = sparseVector(y[wch], wch, length = length(y))

system.time({
  out_s <- moltenNMF:::mNMF_svb(Y, X = X, L = 2,
                                n_epochs = 100, 
                                n_batches = 2000,
                                lr_param = c(20,0.9), 
                                lr_type = "exponential",
                                display_progress = TRUE)
})
#   user  system elapsed 
# 46.812   0.120  47.079 
head.matrix(out_s$shape)
head.matrix(out_s$rate)

V <- out_s$shape/out_s$rate
f <- moltenNMF::product_m(X, V)

plot(as.numeric(dat$Y), f_d,  pch=1, col=rgb(1,0.5,0,0.1))
points(as.numeric(dat$Y), f,  pch=2, col=rgb(0,0.5,1,0.1))
abline(0, 1, col="grey", lty=2)

plot(as.numeric(log1p(dat$Y)), log1p(f_d),  pch=1, col=rgb(1,0.5,0,0.1))
points(as.numeric(log1p(dat$Y)), log1p(f),  pch=2, col=rgb(0,0.5,1,0.1))
abline(0, 1, col="grey", lty=2)

####
#check param
logV <- digamma(out_d$shape)-log(out_d$rate)
What_d <- logV[1:110,2:1]
What_d <- sweep(What_d,1,rowMeans(What_d))

plot(dat$trueW, What_d, col=rgb(1,0,0,0.2))
abline(0,1,lty=2)

Hhat_d <- logV[111:210,2:1]
Hhat_d <- sweep(Hhat_d,1,rowMeans(Hhat_d))
points(dat$trueH,t(Hhat_d), col=rgb(1,0,0,0.2))
abline(0,1,lty=2)

