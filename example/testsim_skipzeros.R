library(ggplot2)
library(parallel)
library(moltenNMF)
library(Matrix)
library(tidyr)
library(dplyr)

L <- 5L
df <- as.data.frame(expand.grid(x1=factor(1:50),
                                x2=factor(1:50),
                                b1=factor(1:2))) %>% 
  mutate(b2= factor(rbinom(n(),1,0.5)))
dim(df)
X <- sparse_onehot(~ ., data=df, binary_dummy = TRUE)
colnames(X)
N <- nrow(X)
D <- ncol(X)
set.seed(2754)
V <- matrix(rgamma(L*D, 0.5, 0.5),D,L)
ord = order(apply(V,2,var), decreasing = TRUE)
V = V[,ord] #reorder by variance
lambda <- product_m.default(X,V)
Y <- rpois(N, lambda)

system.time({
  out_d <- moltenNMF::mNMF_vb.default(Y, X = X,L = L,iter=100,
                                      a=1, b=1,
                                      display_progress = TRUE)
})

wch = which(Y>0)
Y1 = Y[wch]
X1 = slice_rows(X, wch)

system.time({
  out_sb = moltenNMF:::mNMF_svb_batch(Y1, X1, L=L, N=nrow(X), iter=100)
})

Y_sp = sparseVector(Y1, wch, length = length(Y)) 
system.time({
  out_sb2 = moltenNMF:::mNMF_svb_batch(Y_sp, X1, L=L, N=nrow(X), iter=100)
})

plot(out_sb$ELBO[-1],type = "l")
points(out_sb2$ELBO[-1],type = "l", col="royalblue")

V_d <- out_d$shape/out_d$rate
f_d <- moltenNMF::product_m(X, V_d)
V_sb <- out_sb$shape/out_sb$rate
f_sb <- moltenNMF::product_m(X, V_sb)
V_sb2 <- out_sb2$shape/out_sb2$rate
f_sb2 <- moltenNMF::product_m(X, V_sb2)

ggplot()+
  geom_point(aes(x=f_d, y=as.matrix(Y)), alpha=0.25, shape=1)+
  geom_point(aes(x=f_sb, y=as.matrix(Y)), alpha=0.25, colour="royalblue", shape=2)+
  geom_point(aes(x=f_sb2, y=as.matrix(Y)), alpha=0.25, colour="cornflowerblue", shape=2)+
  geom_abline(intercept = 0, slope=1, linetype=2, colour="lightgrey")+
  theme_bw()

#######
#SVB
#######
wch = which(Y>0)
Y1 = Y[wch]
X1 = slice_rows(X, wch)
dim(X1)
system.time({
  out_s <- moltenNMF:::mNMF_svb(Y1, X = X1,
                                N = nrow(X), L = L,
                                n_epochs = 200,
                                n_batches = 100,
                                lr_param = c(15,0.9),
                                lr_type = "exponential",
                                display_progress = TRUE)
})

system.time({
  out_s <- moltenNMF:::mNMF_svb(Y_sp, X = X1,
                                N = nrow(X), L = L,
                                n_epochs = 200,
                                n_batches = 100,
                                lr_param = c(15,0.9),
                                lr_type = "exponential",
                                display_progress = TRUE)
})

plot(out_s$ELBO[-1], type="l")

V_s <- out_s$shape/out_s$rate
f_s <- moltenNMF::product_m(X, V_s)

plot(f_d, as.matrix(Y), pch=1, col=rgb(0,0,0,0.5), cex=0.5, xlab = "fitted")
points(f_s, as.matrix(Y),  pch=1, col=rgb(0,0.5,1,0.5), cex=0.5)
points(as.matrix(Y), lambda,  pch=1, col=rgb(1,0.5,0,0.5), cex=0.5)
abline(0, 1, col="grey", lty=2)

cor(V,V_s)
cor(V,V_d)


rearrange_winner_ord <- function(V,V_s){
  cmat = cor(V, V_s)
  ord =integer(ncol(V))
  ord[1] = which.max(cmat[1,])
  for(i in 2:L){
    ord[i] = which(cmat[i,]==max(cmat[i,-ord[1:(i-1)]]))  
  }
  list(V=V_s[,ord], cor=diag(cmat[,ord]))
}

reV_s = rearrange_winner_ord(V,V_s)
reV_d = rearrange_winner_ord(V,V_d)
# plot(V, reV_d$V)
# points(V, reV_s$V, col="royalblue", pch=2)
plot(V, reV_s$V, col="royalblue", pch=2)

