library(ggplot2)
library(parallel)
library(moltenNMF)
library(Matrix)
library(tidyr)
library(dplyr)

# 3-way tensor
rearrange_winner_ord <- function(V, V_s){
  cmat = cor(V, V_s)
  ord =integer(ncol(V))
  ord[1] = which.max(cmat[1,])
  for(i in 2:L){
    ord[i] = which(cmat[i,]==max(cmat[i,-ord[1:(i-1)]]))  
  }
  list(V=V_s[,ord], cor=diag(cmat[,ord]))
}

simfunc <- function(seed, Xmat, V, L, lambda, settings){
  set.seed(seed)
  Y <- rpois(N, lambda)
  out_d <- moltenNMF::mNMF_vb.default(Y, Xmat, L = L, a = 1, b=1, iter=1000,
                                      display_progress = FALSE)
  V_d <- out_d$shape/out_d$rate
  resV_d = rearrange_winner_ord(V, V_d)
  wch = which(Y>0)
  Y1 = Y[wch]
  X1 = slice_rows(X, wch)
  
  Vcor = vector("list", nrow(settings))
  for(set in seq_len(nrow(settings))){
    lr_param = c(settings$delay[set],
                 settings$forgetting[set])
    out_s <- moltenNMF:::mNMF_svb(Y1, X = X1, N=nrow(X),
                                  L = L,
                                  n_batches  = settings$n_batches[set],
                                  n_epochs = 200,
                                  lr_param = lr_param,
                                  lr_type = "exponential",
                                  display_progress = FALSE)
    V_s = out_s$shape/out_s$rate
    resV_s = rearrange_winner_ord(V_s,V)
    
    Vcor[[set]] = resV_s$cor
  }
  
  list(svb=Vcor, bvb=resV_d$cor)
}

settings = expand.grid(forgetting = c(0.7,0.8,0.9),
                       delay = c(1.5,5,15),
                       n_batches = c(500,1000,2000))

ncols = c(100, 500, 1000)
L <- 5L
df <- as.data.frame(expand.grid(row=factor(1:10),
                                col=factor(1:ncols[1]),
                                depth=factor(1:10)))
dim(df)
X <- sparse_onehot(~ ., data=df)

N <- nrow(X)
D <- ncol(X)
set.seed(1111);V <- matrix(rgamma(L*D, 0.5, 0.5),D,L)
ord = order(apply(V, 2, var), decreasing = TRUE)
V = V[,ord]
lambda <- product_m.default(X, V)

res = simfunc(1, Xmat=X, V=V, L=L, lambda=lambda, settings = settings)

plot(res$svb[[1]], res$bvb)
abline(0,1,lty=2)
