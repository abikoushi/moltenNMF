library(ggplot2)
library(parallel)
library(moltenNMF)
library(Matrix)
library(tidyr)
L <- 5L
# 2-way

settings = expand.grid(forgetting = c(0.9,1),
                       delay = c(1,5,15),
                       n_batches = c(500,1000,2000))

settings[1,, drop=FALSE]$forgetting
#[1] 18  3

set_data_mf <- function(L, nrow, ncol, mu=0){
  W <- matrix(rnorm(nrow*L,0,1), ncol=L)
  H <- matrix(rnorm(L*ncol,0,1), nrow=L)
  W <- sweep(W,1,rowMeans(W)-mu)
  H <- sweep(H,2,colMeans(H)-mu)
  list(lambda=exp(W)%*%exp(H), trueW=exp(W), trueH=exp(H))
}

param <- set_data_mf(L,10,10)

lambda <- param$lambda
simfunc <- function(seed, Xmat, V, L, lambda){
  set.seed(seed)
  Y <- matrix(rpois(length(lambda),lambda),
              nrow(lambda), ncol(lambda))
  Y <- as(Y, "TsparseMatrix")

  out_d <- moltenNMF:::NMF2D_vb(Y, rank = L, 
                                prior_shape = 1, prior_rate = 1, iter=1000,
                                display_progress = FALSE)
  V_d = moltenNMF:::meanV_array(out_d)
  resV_d = rearrange_winner_ord(rbind(V[[1]],V[[2]]), rbind(param$trueW,t(param$trueH)))
  

  out_s <- moltenNMF:::NMF2D_svb(Y, rank = L,
                                n_epochs = 200,
                                n_baches  = 2000,
                                lr_param = c(15,0.9),
                                lr_type = "exponential",
                                display_progress = FALSE)
  
  V_s = moltenNMF:::meanV_array(out_s)
  resV_s = rearrange_winner_ord(rbind(V_s[[1]],V_s[[2]]),
                                rbind(param$trueW,t(param$trueH)))

  return(cbind(resV_d$cor, resV_s$cor))
}



###

cmat = cor(V, V_s)
ord =integer(ncol(V))
ord[1] = which.max(cmat[1,])
for(i in 2:L){
  ord[i] = which(cmat[i,]==max(cmat[i,-ord[1:(i-1)]]))  
}
list(V=V_s[,ord], cor=diag(cmat[,ord]))


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

simfunc <- function(seed, Xmat, V, L, lambda){
  set.seed(seed)
  Y <- rpois(N, lambda)
  out_d <- moltenNMF::mNMF_vb.default(Y, Xmat, L = L, a = 1, b=1, iter=1000,
                                      display_progress = FALSE)
  V_d <- out_d$shape/out_d$rate
  resV_d = rearrange_winner_ord(V, V_d)
  wch = which(Y>0)
  Y1 = Y[wch]
  X1 = slice_rows(X, wch)
  out_s <- moltenNMF:::mNMF_svb(Y1, X = X1,
                                N = nrow(X), L = L,
                                n_epochs = 200,
                                n_batches = 2000,
                                lr_param = c(15,0.9),
                                lr_type = "exponential",
                                display_progress = FALSE)
  
  V_s <- out_s$shape/out_s$rate
  resV_s = rearrange_winner_ord(V,V_d)
  return(cbind(resV_d$cor, resV_s$cor))
}


n_batches = 2000
lr_param = c(15,0.9)


L <- 5L
df <- as.data.frame(expand.grid(row=factor(1:20),
                                col=factor(1:20),
                                depth=factor(1:10)))
dim(df)
X <- sparse_onehot(~ ., data=df)

N <- nrow(X)
D <- ncol(X)
set.seed(1111);V <- matrix(rgamma(L*D, 0.5, 0.5),D,L)
ord = order(apply(V, 2, var), decreasing = TRUE)
V = V[,ord]
lambda <- product_m.default(X, V)

res = lapply(1:3, simfunc, Xmat=X, V=V, L=L, lambda=lambda)

plot(do.call("rbind", res))
abline(0,1,lty=2)
