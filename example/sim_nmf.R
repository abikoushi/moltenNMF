library(moltenNMF)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS)

center0 <- function(X,B){
  i <- attr(X,"assign")
  do.call("rbind", lapply(1:3, function(j)sweep(B[i==j,],2,colMeans(B[i==j,]))))  
}


set_data <- function(L, nrow, ncol){
  df <- expand.grid(row=1:nrow, col=1:ncol, dep=1:10)
  X <- sparse_onehot(~row+col+dep, data=df)
  B <- matrix(rnorm((nrow+ncol+10)*L,0,1),ncol=L)
  B <- center0(X,B)
  y <- rpois(nrow(X), rowSums(exp(X%*%B)))
  list(X=X, y=y, df=df, trueparam =B)
}
set.seed(555)
dat <- set_data(3, 1000, 1000)

X1 <- with(dat, X[y>0,])
y1 <- with(dat, y[y>0])

#attributes(X1) <- attributes(dat$X)

(m <- nrow(dat$X)-nrow(X1))
p <- with(dat, colMeans(X[y==0,]))
m
m/nrow(dat$X)

system.time({
  out1 <- moltenNMF::mNMF_vb(dat$y, dat$X,
                              L = 3,
                              iter = 200)
})
#   user  system elapsed 
#326.887  17.117 344.981 
   
system.time({
  out2 <- moltenNMF:::mNMF_vb_sp(y1, X1,
                                 L = 3,
                                 indices = attr(dat$X,"indices"),
                                 N0 = m,
                                 probx = p,
                                 iter = 200)
})
#   user  system elapsed 
#316.385  16.301 333.673 

f <- moltenNMF::product_m(dat$X, out2$shape/out2$rate)
plot(f, dat$y, pch=".")
abline(0, 1, col="royalblue")

ind <- apply(cor(dat$trueparam,digamma(out2$shape)-log(out2$rate)),1,which.max)
dim(Bhat)
Bhat <- (digamma(out2$shape)-log(out2$rate))[,ind]
Bhat <- center0(dat$X,Bhat)

plot(dat$trueparam,Bhat,pch=".")
abline(0,1,col="royalblue")

plot(out2$ELBO[-1], type = "l")

########
dfest <- data.frame(digamma(out2$shape)-log(out2$rate),
                 axis = attr(dat$X,"assign")) %>% 
  pivot_longer(1:2, values_to = "est") %>% 
  group_by(axis, name) %>% 
  #mutate(est = est - mean(est)) %>% 
  ungroup() %>% 
  mutate(true = c(t(dat$trueparam)))

  
dfest

ggplot(dfest, aes(x=est, y=true))+
  geom_point(alpha=0.3)+
  geom_abline(intercept=0, slope=1, linetype=2)+
  facet_grid(axis ~ name)+
  theme_bw()+labs(x="log-estimates", colour="method")
ggsave("comp.png")
