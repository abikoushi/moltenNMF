library(moltenNMF)
library(parallel)
library(ggplot2)

simmodelselect <- function(seed,X,L){
  set.seed(seed)
  V <- matrix(rgamma(L*ncol(X), 1, 0.01),ncol(X),L)
  lambda <- product_m.default(X,V)
  Y <- rpois(nrow(X), lambda)
  ELBO <- numeric(9)
  for(i in 1:9){
    out <- mNMF_vb.default(Y, X, L = i+1L, a = 0.5, b=1, iter=20000)
    ELBO[i] <- out$ELBO[length(out$ELBO)]
  }
  return(ELBO)
}

df <- as.data.frame(expand.grid(row=factor(1:10),
                                col=factor(1:10),
                                depth=factor(1:10)))
X <- sparse_model_matrix_b(~ . -1, data=df)
simELBOout <- mclapply(1:50, function(i)simmodelselect(i, X=X, L=5),
                       mc.cores = detectCores())
simELBOmat <- simplify2array(simELBOout)
#saveRDS(simELBOmat, file = "./example/simELBOmat5.rds")
table(apply(simELBOmat, 2, which.max)+1)

df <- data.frame(L = rep(2:10, dim(simELBOmat)[2]),
           id = rep(1:dim(simELBOmat)[2], each=9),
           value = c(simELBOmat))

ggplot(df, aes(x=L, y = -log(-value), group=id))+
  geom_line(alpha=0.1) +
  scale_x_continuous(breaks = seq(2, 10, by=3)) +
  theme_classic(18)
ggsave("simELBO5.pdf")
# ggplot(df, aes(x=L, y=value, group=L))+
#   geom_violin(trim = FALSE)+
#   stat_summary(fun.data = mean_sdl, fun.args=list(mult=1))+
#   scale_x_continuous(n.breaks = 9)+
#   theme_classic(18)

#jpeg("~/Desktop/simELBO.jpeg")
barplot(table(apply(simELBOmat, 2, which.max)+1),xlab = "L",ylab = "count",cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5)
#dev.off()
# dfELBO <- as.data.frame(simELBOoutmat) %>% 
#   set_names(2:9) %>% 
#   mutate(id=row_number()) %>% 
#   gather(L,ELBO,-id)

# ggplot(dfELBO,aes(x=L,y=ELBO))+
#   geom_boxplot()+
#   theme_bw(24)