library(NMF)
library(Matrix)
library(VBsNMF)

simNMF <- function(n_row, n_col, den, L){
  maxit <- 100
  meth <- c("brunet", "lee")
  Y <- rsparsematrix(n_row, n_col, den,
                     rand.x= function(n){rpois(n,5)+1},
                     repr = "T")
  Y2 <- as.matrix(Y)
  i0 <- which(colSums(Y2) > 0)
  Y <- Y[,i0]
  Y2 <- Y2[,i0]
  t0 <- system.time({out0 = vb_nmf_pois(Y, rank=L, iter=maxit)})
  t1 <- system.time({out1 = nmf(Y2, rank=L, maxIter=maxit, eps=0,
                                method = meth[1])})
  t2 <- system.time({out2 = nmf(Y2, rank=L, maxIter=maxit, eps=0,
                                method = meth[2])})
  
  fit0 <- basemean(out0)%*%coefmean(out0)
  fit1 <- basis(out1)%*%coef(out1)
  fit2 <- basis(out2)%*%coef(out2)
  
  rmse0 <- sqrt(mean((Y2-fit0)^2))
  rmse1 <- sqrt(mean((Y2-fit1)^2))
  rmse2 <- sqrt(mean((Y2-fit2)^2))
  
  poisloss0 <- mean(dpois(Y2, fit0, log = TRUE))
  poisloss1 <- mean(dpois(Y2, fit1, log = TRUE))
  poisloss2 <- mean(dpois(Y2, fit2, log = TRUE))
  
  res <- matrix(c(t0[3],t1[3],t2[3],
                  rmse0,rmse1,rmse2,
                  poisloss0,poisloss1,poisloss2),3,3)
  row.names(res) <- c("prop","brunet", "lee")
  return(res)  
}

dens <- seq(0.1, 0.5, by=0.1) #proportion of non-zero entries 
n_rows <- c(100, 1000, 5000)
n_cols <- c(2000)
L_rank <- c(2L, 5L, 10L)
repl <- 10 
df <- expand.grid(n_rows=n_rows,n_cols=n_cols, dens=dens, L_rank=L_rank)

results <- vector("list", nrow(df))
pb <- txtProgressBar(0, nrow(df), style = 3)
for(i in 1:nrow(df)){
  set.seed(i)
  A <- array(0, dim=c(3,3,repl))
  for(j in 1:repl){
    A[,,j] <- with(df, simNMF(n_rows[i], n_cols[i], dens[i], L_rank[i]))
  }
  results[[i]] <- A
  setTxtProgressBar(pb,i)
}

A <- simplify2array(results)
saveRDS(A, file = "vsNMF2.rds")

