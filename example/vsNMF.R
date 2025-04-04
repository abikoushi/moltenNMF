library(rliger)
library(NMF)
library(NNLM)
library(Matrix)
library(moltenNMF)
library(dplyr)

zeroprob_norm <- function(mean, sd=1){
  pnorm(0, mean, sd)
}

x <- seq(-10,10, length.out=100)
y <- seq(-10,10, length.out=200)

basefun <- function(x){
  cbind(tanh(x),sin(x),cos(x),sign(x))  
}

softplus <- function(x){log1p(exp(x))}
relu <- function(x){
  x[x<0] <- 0
  x
}

row_norm <- function(V){
  sweep(V, 1, rowSums(V), "/")  
}

B1 <- matrix(rnorm(4*4,0,1),4,4)
B2 <- matrix(rnorm(4*4,0,1),4,4)
V1 <- row_norm(relu(basefun(x)%*%B1))
V2 <- relu(B2%*%t(basefun(y)))
muZ = V1%*%V2
sigma=1
mean(zeroprob_norm(muZ,sigma))

image.default(muZ)

Z = relu(matrix(rnorm(length(x)*length(y),muZ, sigma), length(x), length(y)))
sp = mean(Z==0)
image(x,y,Z)

Z = as(Z, "TsparseMatrix")
image(Z)

system.time2 <- function(expr){
  mem_before <- pryr::mem_used()
  t0 = system.time(expr = expr)
  mem_after <- pryr::mem_used()
  cname = attr(t0, "names")
  t0 = as.data.frame(t(as.numeric(t0)))
  t0 = setNames(t0, cname)
  t0$mem_used = mem_after - mem_before
  return(t0)
}

fit_liger <- function(liger_obj){
  W = getMatrix(liger_obj, slot = "W")
  V = getMatrix(liger_obj, slot = "V")
  H = getMatrix(liger_obj, slot = "H")
  as.matrix((W + V[[1]])%*%H[[1]])
}

extract_row_liger <- function(liger_obj){
  W = getMatrix(liger_obj, slot = "W")
  V = getMatrix(liger_obj, slot = "V")
  as.matrix(W + V[[1]])
}

fit_mnmf <- function(out){
  Vhat_t <- moltenNMF:::meanV_array(out)
  Vhat_t[[1]]%*%t(Vhat_t[[2]])
}

xlogy <- function(x,y){
  ifelse(x==0 & y==0, 0, x*log(y))  
}

poisloss <- function(obs, fit){
  -mean(xlogy(obs,fit)-fit-lgamma(obs+1))
}

rmse <- function(obs, fit){
  sqrt(mean((as.matrix(obs) - fit)^2))
}

maxit=10
rank = 4
t0 <- system.time2({
  decomp0 <- moltenNMF:::NMF2D_vb(Z, rank=rank, iter=maxit,
                                               display_progress = FALSE)
  })
t0
t1 <- system.time2({
  decomp1 <- NMF::nmf(as.matrix(Z), rank = rank, maxIter=maxit, eps=0)
})
t1

t2 = system.time2({
  decomp2 <- nnmf(as.matrix(Z), k = rank, rel.tol = 0, max.iter = maxit, loss = "mkl", verbose=0, method = "lee")
})
t2

t3 <- system.time2({
  rownames(Z) <- paste0("gene", 1:nrow(Z))
  colnames(Z) <- paste0("cell", 1:ncol(Z))
  liger_obj = createLiger( list(sample1=Z) )
  liger_obj@datasets$sample1@normData <- liger_obj@datasets$sample1@rawData
  liger_obj@datasets$sample1@scaleData <- liger_obj@datasets$sample1@normData
  liger_obj@varFeatures <- rownames(liger_obj@datasets[[1]])
  liger_obj <- runINMF(liger_obj,
                      k=rank, nCores = 1,
                       nIteration = maxit, 
                       verbose = FALSE)
})

Vhat_t <- moltenNMF:::meanV_array(decomp0)
fit_0 <- Vhat_t[[1]]%*%t(Vhat_t[[2]])
fit_1 <- NMF::basis(decomp1)%*%NMF::coef(decomp1)
fit_2 <- decomp2$W%*%decomp2$H
fit_3 <- fit_liger(liger_obj)

Vhat0 = row_norm(Vhat_t[[1]])
Vhat1 = row_norm(NMF::basis(decomp1))
Vhat2 = row_norm(decomp2$W)
Vhat3 = row_norm(extract_row_liger(liger_obj))
matplot(V1, type="l",col="black", ylim=c(0,1))
matlines(Vhat0, type="l",col="gold")
matlines(Vhat1, type="l",col="cornflowerblue")
matlines(Vhat2, type="l",col="orangered")
matlines(Vhat3, type="l",col="darkviolet")


t0 = mutate(t0, poisloss = poisloss(as.matrix(Z), fit_0), 
             rmse = rmse(as.matrix(Z), fit_0), method="proposed")
t1 = mutate(t1, poisloss = poisloss(as.matrix(Z), fit_1),
             rmse = rmse(as.matrix(Z), fit_1), method="NMF")
t2 = mutate(t2, poisloss = poisloss(as.matrix(Z), fit_2),
             rmse = rmse(as.matrix(Z), fit_2) , method="NNLM")
t3 = mutate(t3, poisloss = poisloss(as.matrix(Z), fit_3),
            rmse = rmse(as.matrix(Z), fit_3), method="rliger")

bm = rbind(t0,t1,t2,t3)
bm

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
#saveRDS(A, file = "vsNMF2.rds")

