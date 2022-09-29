library(ggplot2)
library(parallel)
library(moltenNMF)
#library(Matrix)
#library(Rcpp)

simfunc <- function(seed,V,L,varind,perm_L){
  set.seed(seed)
  Y <- rpois(N, lambda)
  out <- mNMF_vb.default(Y, X, L = L, a = 0.5, b=1, iter=10000)
  Vhat <- out$shape/out$rate
  #3 is number of variables
  cmat <- sapply(1:3,function(j){
    apply(perm_L,1,function(ivec)cor(log(as.vector(Vhat[varind==j,ivec])),
                                   log(as.vector(V[varind==j,]))))
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
X <- sparse_model_matrix_b(~ . -1, data=df)
N <- nrow(X)
D <- ncol(X)
set.seed(1111);V <- matrix(rgamma(L*D, 1, 0.01),D,L)
lambda <- product_m.default(X,V)
#3 is number of variables
varind <- rep(1:3,diff(attr(X, "indices")))
## plot(out$ELBO, type = "l")
system.time({
  simout <- mclapply(1:100, function(i)simfunc(i,V,L,varind,perm_L), mc.cores = detectCores())  
})

save(simout,V,X,file = "./example/sim_poisson.Rdata")

Vmean <- apply(simplify2array(lapply(simout,function(x)x$Vhat)),1:2,mean)
Vsd <- apply(simplify2array(lapply(simout,function(x)x$Vhat)),1:2,sd)
length(simout)
dfV <- data.frame(true=as.vector(V),
                  mean=as.vector(Vmean),
                  se=as.vector(Vsd),
                  component=rep(1:L, each=nrow(V)),
                  variable=rep(varind,L))

p1 <- ggplot(dfV,aes(x=true,y=mean,ymin=mean-se,ymax=mean+se))+
  stat_smooth(method = lm,formula=y ~ x - 1, linetype=2, se=FALSE, colour="lightgrey")+
  geom_pointrange(shape=1)+
  scale_x_continuous(n.breaks = 4)+
  facet_wrap(variable~component, scales = "free", labeller = label_both, ncol=5)+
  xlab("true value")+ylab("esitmates")+
  theme_classic(14) +theme(strip.background = element_rect(colour = NA, fill="grey95"))
print(p1)
#ggsave(p1, filename = "simV.pdf")

