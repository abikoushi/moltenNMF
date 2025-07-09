library(ggplot2)
library(parallel)
library(moltenNMF)
library(Matrix)
library(tidyr)
library(dplyr)
library(foreach)
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

set.seed(1111);V <- matrix(rgamma(L*D, 0.5, 0.5),D,L)
ord = order(apply(V, 2, var), decreasing = TRUE)
V = V[,ord]
lambda <- product_m.default(X, V)
i=1
set.seed(i)
Y <- rpois(N, lambda)
out_d <- moltenNMF::mNMF_vb.default(Y, X, L = L, a = 1, b=1, iter=100,
                                    display_progress = FALSE)
V_d <- out_d$shape/out_d$rate
resV_d = rearrange_winner_ord(V, V_d)
wch = which(Y>0)
Y1 = Y[wch]
X1 = slice_rows(X, wch)
ncols = c(100, 500, 1000)
settings = expand.grid(forgetting = c(0.7,0.8,0.9),
                       delay = c(1.5,5,15),
                       n_batches = c(500,1000,2000),
                       rep = 1:10)

res = foreach(i = 1, 
         .export = c("Y1","X1","X","rearrange_winner_ord","V","settings")) %do%{
  L <- 5L
  df <- as.data.frame(expand.grid(row=factor(1:50),
                                  col=factor(1:ncols[3]),
                                  depth=factor(1:2)))
  
  X <- sparse_onehot(~ ., data=df)
  N <- nrow(X)
  D <- ncol(X)
  lr_param = c(settings$delay[i], settings$forgetting[i])
  out_s <- moltenNMF:::mNMF_svb(Y1, X = X1, N=nrow(X),
                                L = L,
                                n_batches  = settings$n_batches[i],
                                n_epochs = 10,
                                lr_param = lr_param,
                                lr_type = "exponential",
                                display_progress = FALSE)
    V_s = out_s$shape/out_s$rate
    resV_s = rearrange_winner_ord(V_s,V)
    resV_s$cor
}


# f_s <- moltenNMF::product_m(X, V_s)

# ggplot(data = NULL)+
#   geom_abline(slope = 1, intercept = 0, colour="lightgrey")+
#   geom_point(aes(x=c(as.matrix(Y)), y=c(f_s)), alpha=0.1, size=1)+
#   geom_bin_2d(aes(x=c(as.matrix(Y)), y=c(lambda),fill=after_stat(log10(count))), alpha=0.25)+
#   scale_fill_viridis_c()+
#   theme_bw(16)
