library(moltenNMF)
library(dplyr)
library(Matrix)
Titanicdf <- mutate(as.data.frame(Titanic), 
                    Class=factor(Class,levels=c("3rd","2nd","1st","Crew")))
Titanicdf2 <- Titanicdf
ds <- dim(Titanicdf[,-5])
ind <- cbind(sample.int(ds[1],5,replace = TRUE),
             sample.int(ds[2],5,replace = TRUE))
Titanicdf2[ind] <- NA
f <- Freq ~ Survived+Class+Sex+Age-1
X <- sparse_model_matrix_b(f, Titanicdf2)
M <- which(is.na(Titanicdf2),arr.ind = TRUE)-1L
out <- mNMF_vb.default(Titanicdf2$Freq, X = X, L=2, pos_missing = M, iter=10)

Titanicdf[ind]
colnames(out$Xprob) <- colnames(X)
out$Xprob
Titanicdf[ind]
