library(parallel)
library(moltenNMF)
library(Matrix)
library(foreach)
#library(doParallel)
#library(snow)
#library(Rmpi)

# cl <- makeCluster(3, type='MPI')
# registerDoParallel(cl)
# 3-way tensor

rearrange_winner_ord <- function(V,V_s){
  cmat = cor(V, V_s)
  ord =integer(ncol(V))
  ord[1] = which.max(cmat[1,])
  for(i in 2:L){
    ord[i] = which(cmat[i,]==max(cmat[i,-ord[1:(i-1)]]))  
  }
  list(V=V_s[,ord], cor=diag(cmat[,ord]))
}

simfunc <- function(seed, V, L, lambda, settings){
  set.seed(seed)
  wch = which(Y>0)
  Y1 = Y[wch]
  X1 = slice_rows(X, wch)
  out_d <- moltenNMF::mNMF_vb.default(Y, X, L = L, a = 1, b=1, 
                                      iter=200,
                                      display_progress = FALSE)
  V_d <- out_d$shape/out_d$rate
  resV_d = rearrange_winner_ord(V, V_d)
  Vcor = vector("list", nrow(settings))
  for(set in seq_len(nrow(settings))){
    lr_param = c(settings$delay[set],
                 settings$forgetting[set])
    out_s <- moltenNMF:::mNMF_svb_batch(Y1, X1, L = L, 
                                        N=length(Y),
                                        n_epochs = 200,
                                        lr_param = lr_param,
                                        lr_type = "exponential",
                                        a = 1, b=1,
                                        display_progress = FALSE)
    V_s <- out_s$shape/out_s$rate
    resV_s = rearrange_winner_ord(V, V_s)
    Vcor[[set]] = resV_s$cor
  }
  list(svb=Vcor, bvb=resV_d$cor)
}

?foreach
res = foreach(i = 1:3,
              .packages = c("Matrix","moltenNMF"),
              .export = c("rearrange_winner_ord", "simfunc")) %dopar% {
  settings = expand.grid(delay = c(0.7,0.8,0.9,1.5),
                         forgetting = c(0.7,0.8,0.9),
                         rep = 1:10)
  L <- 5L
  ncols = c(100, 500, 1000)
  df <- as.data.frame(expand.grid(row=factor(1:50),
                                  col=factor(1:ncols[i]),
                                  depth=factor(1:2)))
  X <- sparse_onehot(~ ., data=df)
  N <- nrow(X)
  D <- ncol(X)
  set.seed(1111);V <- matrix(rgamma(L*D, 0.5, 0.5),D,L)
  ord = order(apply(V, 2, var), decreasing = TRUE)
  V = V[,ord]
  lambda <- product_m.default(X, V)  
  Y <- rpois(N, lambda)
  ressim = simfunc(i, V, L, lambda, settings)
  saveRDS(ressim, paste0("simmnmf_",ncols[i],".rds"))
}


# V_s <- out_s$shape/out_s$rate
# resV_s = rearrange_winner_ord(V, V_s)
# f_s <- moltenNMF::product_m(X, V_s)
# plot(out_s$ELBO[-1], type = "l")
# ggplot(data = NULL)+
#   geom_abline(slope = 1, intercept = 0, colour="lightgrey")+
#   geom_point(aes(x=c(as.matrix(Y)), y=c(f_s)), alpha=0.1, size=1)+
#   geom_bin_2d(aes(x=c(as.matrix(Y)), y=c(lambda),fill=after_stat(log10(count))), alpha=0.25)+
#   scale_fill_viridis_c()+
#   theme_bw(16)
