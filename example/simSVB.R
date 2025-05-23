library(ggplot2)
library(parallel)
library(moltenNMF)
library(Matrix)
library(tidyr)
L <- 5L
simfunc <- function(seed, Xmat, V, L, lambda){
  set.seed(seed)
  Y <- rpois(N, lambda)
  out_d <- moltenNMF::mNMF_vb.default(Y, Xmat, L = L, a = 1, b=1, iter=1000)
  V_d <- out_d$shape/out_d$rate
  resV_d = rearrange_winner_ord(V,V_d)
  out_s <- moltenNMF:::mNMF_svb(Y1, X = X1,
                                N = nrow(X), L = L,
                                n_epochs = 200,
                                n_batches = 2000,
                                lr_param = c(15,0.9),
                                lr_type = "exponential",
                                display_progress = FALSE)
  
  V_s <- out_s$shape/out_s$rate
  resV_s = rearrange_winner_ord(V,V_d)
  return(resV_d$cor)
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
ord = order(apply(V,2,var), decreasing = TRUE)
V = V[,ord]
lambda <- product_m.default(X,V)

res = lapply(1:3, simfunc, Xmat=X, V=V, L=L, lambda=lambda)
res

#condition 
