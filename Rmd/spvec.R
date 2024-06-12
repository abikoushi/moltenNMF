library(moltenNMF)
library(Matrix)
library(dplyr)
library(ggplot2)
set.seed(2222)
lambda<-rexp(60,1)
df <- expand.grid(row=1:20,col=1:20,dep=1:20)
X <- sparse_onehot(~row+col+dep, data=df)
y <- rpois(8000,exp(as.numeric(X%*%log(lambda))))
hist(y)

X1 <- X[y>0,]
y1 <- y[y>0]

(m <- nrow(X)-nrow(X1))


system.time({
  out1 <- moltenNMF::mNMF_vb.default(y,X,1,iter = 500)
})
dim(X)
m
p <- colMeans(X[y==0,])
#p <- (colSums(X[y==0,])+1)/(m+2) #Laplace smoothing
system.time({
  out2 <- moltenNMF:::mNMF_vb_sp(y1,X1,1,
                                 indices = attr(X,"indices"),
                                 N0 = m,
                                 probx = p,
                                 iter = 500)
})


plot(out1$ELBO[-1], type = "l")

plot(out2$ELBO[-1], type = "l")

out2$ELBO

mean((digamma(out1$shape)-log(out1$rate))^2)
mean((digamma(out2$shape)-log(out2$rate))^2)

df <- data.frame(true=lambda,
                 est1=out1$shape/out1$rate,
                 est2=out2$shape/out2$rate,
                 logest1=digamma(out1$shape)-log(out1$rate),
                 logest2=digamma(out2$shape)-log(out2$rate)) %>% 
  mutate(win = if_else((log(true)-logest1)^2 < (log(true)-logest2)^2,"conventional","proposed"))

ggplot(df,aes(x=true))+
  #geom_linerange(aes(ymin=est1,ymax=est2, colour=win),alpha=0.7)+
  geom_point(aes(y=est1,colour="conventional"),alpha=0.7)+
  geom_point(aes(y=est2,colour="proposed"),alpha=0.7)+
  geom_abline(intercept=0,slope=1,linetype=2)+
  theme_bw()+labs(y="estimates",colour="method")

ggplot(df,aes(x=log(true)))+
  #geom_linerange(aes(ymin=logest1,ymax=logest2, colour=win),alpha=0.7)+
  geom_point(aes(y=logest1,colour="conventional"),alpha=0.7)+
  geom_point(aes(y=logest2,colour="proposed"),alpha=0.7)+
  geom_abline(intercept=0,slope=1,linetype=2)+
  theme_bw()+labs(y="estimates",colour="method")

#######
titanicdf <- as.data.frame(Titanic) %>% 
  dplyr::filter(Survived=="Yes")

y <- titanicdf$Freq
X <- sparse_onehot(~Class+Sex+Age, data=titanicdf)
system.time({
  out1 <- moltenNMF::mNMF_vb.default(y,X,2,iter = 200)
})
y0 <- titanicdf$Freq[titanicdf$Freq>0L]
X0 <- sparse_onehot(~Class+Sex+Age, data=titanicdf[titanicdf$Freq>0L,])

N0 <- sum(titanicdf$Freq==0L)
a <- attr(X,"assign")
tab <- colSums(X[titanicdf$Freq==0L,])

#probx <- unlist(tapply(tab, a, function(x)x/sum(x)))

system.time({
  out2 <- moltenNMF:::mNMF_vb_sp(y0, X0, 2,
                                 N0 = N0,
                                 probx = tab/N0,
                                 b=0.01,
                                 iter = 200)
})
v1 <- out1$shape/out1$rate
v1 <- v1/rowSums(v1)
v2 <- out2$shape/out2$rate
v2 <- v2/rowSums(v2)
pheatmap::pheatmap(v1)
pheatmap::pheatmap(v2)
plot(out2$ELBO[-1],type="l")

plot(v1,v2)
abline(0,1,lty=2)

