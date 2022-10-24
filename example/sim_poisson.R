library(ggplot2)
library(parallel)
library(moltenNMF)
library(Matrix)
library(tidyr)
#library(Rcpp)

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
df <- as.data.frame(expand.grid(row=factor(1:10),
                                col=factor(1:10),
                                depth=factor(1:10)))
#X <- sparse_model_matrix_b(~ . -1, data=df)
# N <- nrow(X)
# D <- ncol(X)
X <- sparse_onehot(~ ., data=df)
colnames(X)
N <- nrow(X)
D <- ncol(X)
set.seed(1111);V <- matrix(rgamma(L*D, 0.8, 0.01),D,L)
set.seed(1111);V <- matrix(rgamma(L*D, 1.5, 0.01),D,L)
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
