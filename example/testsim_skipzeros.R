library(ggplot2)
library(parallel)
library(moltenNMF)
library(Matrix)
library(tidyr)
#library(Rcpp)
L <- 5L
perm_L <- gtools::permutations(L,L)
df <- as.data.frame(expand.grid(x=factor(1:5),
                                row=factor(1:10),
                                col=factor(1:10),
                                depth=factor(1:10)))
                                
dim(df)
#X <- sparse_model_matrix_b(~ . -1, data=df)
X <- sparse_onehot(~ ., data=df)
# colnames(X)
N <- nrow(X)
D <- ncol(X)
set.seed(1645)
V <- matrix(rgamma(L*D, 0.5, 0.5),D,L)
ord = order(apply(V,2,var), decreasing = TRUE)
V = V[,ord]
lambda <- product_m.default(X,V)
Y <- rpois(N, lambda)

system.time({
  out_d <- moltenNMF::mNMF_vb.default(Y, X = X,L = L,iter=100,
                                   display_progress = TRUE)
})
V_d <- out_d$shape/out_d$rate
f_d <- moltenNMF::product_m(X, V_d)

wch = which(Y>0)
Y1 = Y[wch]
X1 = slice_rows(X, wch)
X0 = X[-wch,]
diff(X0@p)
system.time({
  out_s <- moltenNMF:::mNMF_vb_sp(Y1, X = X1,
                                  xp0 = X0@p, L = L,
                                  iter=100,
                                  display_progress = TRUE)
})
V_s <- out_s$shape/out_s$rate
f_s <- moltenNMF::product_m(X, V_s)

plot(out_s$ELBO[-1],type="l")

 # plot(f_s,  as.matrix(Y), pch=1, col=rgb(0,0,0,0.2), cex=1)
 # abline(0, 1, col="grey", lty=2)
plot(f_d, as.matrix(Y),  pch=1, col=rgb(0,0,0,0.3), cex=1, xlab="fitted")
points(f_s,  as.matrix(Y), pch=1, col=rgb(0,0.5,1,0.3), cex=1)
# points(as.matrix(Y), lambda,  pch=1, col=rgb(1,0.5,0,0.5), cex=0.5)
abline(0, 1, col="grey", lty=2)

# pop = c(rep(1,10), rep(0,5))
# wch = numeric(10000)
# for(i in 1:10000){
#   s = sample(pop)
#   wch[i] = which.max(s==1)  
# }
# 
# tab = table(wch)/10000
# plot(tab)
# points(1:10, dgeom(1:10-1, 10/15), type = "b")


#######
#SVB
#######
# probX0 = colMeans(X[-wch,])
# N0 = nrow(X)-length(wch)
# system.time({
#   out_s <- moltenNMF:::mNMF_svb_sp(Y1, X = X1,
#                                    N = nrow(X), probX0 = probX0, L = L,
#                                    n_epochs = 100, 
#                                    n_batches = 100,
#                                    lr_param = c(19,0.9), 
#                                    lr_type = "exponential",
#                                    display_progress = FALSE)
# })
# V_s <- out_s$shape/out_s$rate
# f_s <- moltenNMF::product_m(X, V_s)

X0 = X[-wch,]
N0 = nrow(X)-length(wch)
system.time({
  out_s <- moltenNMF:::mNMF_svb_sp(Y1, X = X1,
                                   N = nrow(X), L = L,
                                   n_epochs = 200,
                                   n_batches = 100,
                                   lr_param = c(15,0.9),
                                   lr_type = "exponential",
                                   display_progress = TRUE)
})
V_s <- out_s$shape/out_s$rate
f_s <- moltenNMF::product_m(X, V_s)

plot(f_d, as.matrix(Y), pch=1, col=rgb(0,0,0,0.5), cex=0.5, xlab = "fitted")
points(f_s, as.matrix(Y),  pch=1, col=rgb(0,0.5,1,0.5), cex=0.5)
points(as.matrix(Y), lambda,  pch=1, col=rgb(1,0.5,0,0.5), cex=0.5)
abline(0, 1, col="grey", lty=2)
plot(out_s$ELBO[-1], type="l")

cor(V,V_s)
cor(V,V_d)

cmat = cor(V, V_s)
ord =integer(L)
ord[1] = which.max(cmat[1,])
for(i in 2:L){
  ord[i] = which(cmat[i,]==max(cmat[i,-ord[1:(i-1)]]))  
}

plot(V, V_s[,ord])
diag(cmat[,ord])

cmat = cor(V, V_d)
ord =integer(L)
ord[1] = which.max(cmat[1,])
for(i in 2:L){
  ord[i] = which(cmat[i,]==max(cmat[i,-ord[1:(i-1)]]))  
}

plot(V, V_d[,ord])
diag(cmat[,ord])


simfunc <- function(seed,V,L,perm_L){
  set.seed(seed)
  Y <- rpois(N, lambda)
  out <- mNMF_vb.default(Y, X, L = L, a = 0.5, b=1, iter=1000)
  Vhat <- out$shape/out$rate
  varind <- attr(X,"indices")
  #3 is number of variables
  cmat <- sapply(1:3,function(j){
    apply(perm_L,1,function(ivec)cor(log(as.vector(Vhat[,ivec][(varind[j]+1):varind[j+1],])),
                                     log(as.vector(V[(varind[j]+1):varind[j+1],]))))
  })
  cv <- rowMeans(cmat)
  wch <- which.max(cv)
  ord <- perm_L[wch,]
  list(Vhat = Vhat[,ord], cv = cv[wch])
}
L <- 5L
perm_L <- gtools::permutations(L,L)
df <- as.data.frame(expand.grid(row=factor(1:20),
                                col=factor(1:20),
                                depth=factor(1:10)))
#X <- sparse_model_matrix_b(~ . -1, data=df)
# N <- nrow(X)
# D <- ncol(X)
X <- sparse_onehot(~ ., data=df)
colnames(X)
N <- nrow(X)
D <- ncol(X)
set.seed(1111);V <- matrix(rgamma(L*D, 0.5, 0.5),D,L)
ord = order(apply(V,2,var), decreasing = TRUE)
V = V[,ord]
lambda <- product_m.default(X,V)
#3 is number of variables
## plot(out$ELBO, type = "l")
system.time({
  simout <- mclapply(1:100, function(i)simfunc(i,V,L,perm_L), mc.cores = detectCores())  
})

#save(simout,V,X,file = "./example/sim_poisson.Rdata")
simout[[1]]
Vmean <- apply(simplify2array(lapply(simout,function(x)x$Vhat)),1:2,mean)
Vsd <- apply(simplify2array(lapply(simout,function(x)x$Vhat)),1:2,sd)
dim(Vmean)
dim(V)
dfV <- data.frame(true = as.vector(V),
                  vbar = as.vector(Vmean),
                  se = as.vector(Vsd),
                  component=rep(1:L, each=nrow(V)),
                  varname = rep(colnames(X), L)) %>% 
  separate(varname,c("variable","id"), sep="_")


p1 <- ggplot(dfV,aes(x=true,y=vbar,ymin=vbar-se,ymax=vbar+se, group=paste(variable,component)))+
  stat_smooth(method = lm, formula=y ~ x - 1, se=FALSE, colour="lightgrey")+
  geom_pointrange(shape=1)+
  scale_x_continuous(n.breaks = 3)+
  facet_wrap(~paste(variable,component), scales = "free", ncol=5)+
  xlab("true value")+ylab("esitmates")+
  theme_classic(14)+theme(strip.background = element_blank())
print(p1)
ggsave(p1, filename = "simV.png")
