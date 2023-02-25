library(dplyr)
library(Matrix)
library(moltenNMF)

dat <- mutate(as.data.frame(Titanic), Class=factor(Class,levels=c("3rd","2nd","1st","Crew")))

f <- Freq ~ Survived+Class+Sex+Age
dat <- model.frame(f, data=dat)
X <- sparse_onehot(f, data=dat)
indices <- attr(X,"indices")
y <- model.response(dat)
X_r <- as(X,"RsparseMatrix")
out <- moltenNMF:::doSVB_pois(y, X_r@j, X_r@p, indices, X_r@Dim[2],
                              batch_size = 10, L=2, n_epoch = 1,
                              a=1,
                              b=1)
uv <- sample.int(X@Dim[1],10)
V <- rexp(32)
moltenNMF:::myprod_r_i(N = X@Dim[1], xj = X_r@j, xp = X_r@p, id = uv, lam = matrix(V))
