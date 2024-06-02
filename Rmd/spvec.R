library(moltenNMF)
library(Matrix)
titanicdf <- as.data.frame(Titanic)
y <- as(titanicdf$Freq,"sparseVector")
X <- sparse_onehot(~Class+Sex+Age, data=titanicdf)
system.time({
  out1 <- moltenNMF::mNMF_vb.default(y,X,2,iter = 100000)  
})
system.time({
  out2 <- moltenNMF::mNMF_vb.default(titanicdf$Freq,X,2,iter = 100000)
})

plot(out1$ELBO[-1], type = "l", col=1)
lines(out2$ELBO[-1], type = "l", col=2)

image(out1$shape/out1$rate)
image(out2$shape/out2$rate)

