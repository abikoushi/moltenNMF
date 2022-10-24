library(moltenNMF)
library(dplyr)
library(Matrix)
f <- Freq ~ . -1
Titanicdf <- mutate(as.data.frame(Titanic), 
                    Class=factor(Class,levels=c("3rd","2nd","1st","Crew")))
X <- sparse_onehot(f, Titanicdf)
dim(X)
Titanicdf2 <- Titanicdf
ds <- dim(Titanicdf[,-5])
ind1 <- sort(sample.int(ds[1],5))
ind2 <- sort(sample.int(ds[2],5,replace = TRUE))
Titanicdf2[cbind(ind1,ind2)] <- NA
M <- which(is.na(Titanicdf2), arr.ind = TRUE)
out <- mNMF_vb.default(Titanicdf2$Freq, X = X, L=2, pos_missing = M, iter=20)
plot(out$Xprob, X[M[,1],])
plot(out$ELBO, type = "l")
