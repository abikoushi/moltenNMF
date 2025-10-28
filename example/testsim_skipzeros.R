library(moltenNMF)
library(Matrix)
library(tidyr)
library(dplyr)
library(ggplot2)

moltenNMF:::check_arma64()

rearrange_winner_ord <- function(V,V_s){
  cmat = cor(V, V_s)
  ord =integer(ncol(V))
  ord[1] = which.max(cmat[1,])
  for(i in 2:L){
    ord[i] = which(cmat[i,]==max(cmat[i,-ord[1:(i-1)]]))  
  }
  list(V=V_s[,ord], cor=diag(cmat[,ord]))
}

i = 1
L <- 5L
ncols = c(100, 500, 1000)
df <- as.data.frame(expand.grid(row=factor(1:100),
                                col=factor(1:ncols[i]),
                                depth=factor(1:2)))

X <- sparse_onehot(~ ., data=df)
N <- nrow(X)
D <- ncol(X)
set.seed(1234);V <- matrix(rlnorm(L*D, 0, 1),D,L)
ord = order(apply(log(V), 2, var), decreasing = TRUE)
V = V[,ord]
lambda <- product_m.default(X, V)  
Y <- rpois(N, lambda)
wch = which(Y>0)

Y1 = Y[wch]
X1 = slice_rows(X, wch)

system.time({
  out_s <- moltenNMF:::mNMF_svb(Y1, X = X1,
                                N = nrow(X), L = L,
                                n_batches = 10000,
                                n_epochs = 200,
                                lr_param = c(15,0.9),
                                lr_type = "exponential",
                                M_max = 100,
                                display_progress = TRUE)
})
head(Y1)



#is.list(out_s$shape/out_s$rate)
plot(out_s$ELBO[-1], type = "l")

f=product_m(X1, out_s$shape/out_s$rate)
plot(f, Y1, pch=".", col=rgb(0,0,0,0.5))
abline(0,1,col="royalblue")
# corVs = rearrange_winner_ord(log(V),
#                              V_s = digamma(out_s$shape) - log(out_s$rate))
corVs = rearrange_winner_ord(V,
                             V_s = out_s$shape/out_s$rate)
corVs$cor
plot(log(V),log(corVs$V))

plot(corVs$V[,5])

####

L <- 4L
df1 <- as.data.frame(expand.grid(x1=factor(1:100),
                                x2=factor(1:100)))

df2 <- as.data.frame(expand.grid(x1=factor(101:201),
                                x2=factor(201:211)))
df = mutate(rbind(df1,df2))

X0 <- sparse_onehot(~ ., data=df)

H=matrix(rbinom(nrow(X0)*10,1,0.2),nrow(X0))

X = append_new(X0,H)

N <- nrow(X)
D <- ncol(X)
set.seed(55)
V <- matrix(rgamma(L*D, 0.5, 0.5), D, L)
ord = order(apply(V,2,var), decreasing = TRUE)
V = V[,ord] #reorder by variance
lambda <- product_m.default(X,V)
Y <- rpois(N, lambda)
#######
#SVB
#######
wch = which(Y>0)
Y1 = Y[wch]
X1 = slice_rows(X, wch)
length(Y1)
range(unique(X1@p))
dim(X1)
#rho = length(Y1)/nrow(X)
#1/rho
#(1-rho)/(rho^2)
#length(Y1)

system.time({
  out_s <- moltenNMF:::mNMF_svb(Y1, X = X1,
                                N = nrow(X), L = L,
                                n_epochs = 100,
                                n_batches = 100,
                                lr_param = c(1,0.9),
                                lr_type = "exponential",
                                M_max = 10,
                                display_progress = TRUE)
})

# system.time({
#   out_s2 <- moltenNMF:::mNMF_svb_batch(Y1, X = X1,
#                                 N = nrow(X), L = L,
#                                 n_epochs = 500,
#                                 lr_param = c(1,0.9),
#                                 lr_type = "exponential",
#                                 display_progress = TRUE)
# })

plot(out_s$ELBO[-1], type="l", xlim = c(0,500))
lines(out_s2$ELBO[-1], type="l", col="royalblue")

V_s <- out_s$shape/out_s$rate
f_s <- moltenNMF::product_m(X, V_s)


V_s2 <- out_s2$shape/out_s2$rate
f_s2 <- moltenNMF::product_m(X, V_s2)

ggplot(data = NULL)+
  geom_abline(slope = 1, intercept = 0, colour="lightgrey")+
  geom_bin2d(aes(x=c(as.matrix(Y)), y=c(lambda), fill = after_stat(log10(count))), alpha = 0.2)+
  geom_point(aes(x=c(as.matrix(Y)), y=c(f_s)), alpha=0.5, size=0.5,colour="orangered")+
  geom_point(aes(x=c(as.matrix(Y)), y=c(f_s2)), alpha=0.5, size=0.5, colour="cornflowerblue")+
  scale_fill_viridis_c()+
  theme_bw(16)

####
#batch
####

system.time({
  out_d <- moltenNMF::mNMF_vb.default(Y, X = X,L = L,iter=500,
                                      a=1, b=1,
                                      display_progress = TRUE)
})


wch = which(Y>0)
Y_sp = sparseVector(Y1, wch, length = length(Y)) 
system.time({
  out_sb = moltenNMF:::mNMF_vb.default(Y_sp, X, L=L, iter=100)
})

plot(out_d$ELBO[-1],type = "l")
plot(out_sb$ELBO[-1],type = "l", col="royalblue")


V_d <- out_d$shape/out_d$rate
f_d <- moltenNMF::product_m(X, V_d)
V_sb <- out_sb$shape/out_sb$rate
f_sb <- moltenNMF::product_m(X, V_sb)

ggplot()+
  geom_point(aes(x=f_s, y=as.matrix(Y)), alpha=0.25, shape=1)+
  geom_point(aes(x=f_sb, y=as.matrix(Y)), alpha=0.25, colour="royalblue", shape=2)+
  geom_abline(intercept = 0, slope=1, linetype=2, colour="lightgrey")+
  theme_bw()

rearrange_winner_ord <- function(V,V_s){
  cmat = cor(V, V_s)
  ord =integer(ncol(V))
  ord[1] = which.max(cmat[1,])
  for(i in 2:L){
    ord[i] = which(cmat[i,]==max(cmat[i,-ord[1:(i-1)]]))  
  }
  list(V=V_s[,ord], cor=diag(cmat[,ord]))
}

ba05 = rgb(0,0,0,0.5)
reV_s = rearrange_winner_ord(V,V_s)
plot(log(V), log(reV_s$V), pch=1, col=ba05)
abline(0,1, col="lightgrey")

reV_d = rearrange_winner_ord(V,V_d)
plot(log(V), log(reV_d$V), col=ba05)
abline(0,1, col="lightgrey")

rownames(V_s)
