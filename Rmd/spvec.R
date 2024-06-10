library(moltenNMF)
library(Matrix)
titanicdf <- as.data.frame(Titanic)

y <- titanicdf$Freq
X <- sparse_onehot(~Class+Sex+Age+Survived, data=titanicdf)
system.time({
  out1 <- moltenNMF::mNMF_vb.default(y,X,2,iter = 100)
})
y0 <- titanicdf$Freq[titanicdf$Freq>0L]
X0 <- sparse_onehot(~Class+Sex+Age+Survived, data=titanicdf[titanicdf$Freq>0L,])

N0 <- sum(titanicdf$Freq==0L)
a <- attr(X,"assign")
tab <- colSums(X[titanicdf$Freq==0L,])

#probx <- unlist(tapply(tab, a, function(x)x/sum(x)))

system.time({
  out2 <- moltenNMF:::mNMF_vb_sp(y0, X0, 2,
                                 N0 = N0,
                                 probx = tab/N0,
                                 #b=0.01,
                                 iter = 200)
})
v1 <- out1$shape/out1$rate
v1 <- v1/rowSums(v1)
v2 <- out2$shape/out2$rate
v2 <- v2/rowSums(v2)
pheatmap::pheatmap((v1))
pheatmap::pheatmap((v2))
plot(out2$ELBO[-1],type="l")

plot(v1,v2[,1:2])
abline(0,1,lty=2)

