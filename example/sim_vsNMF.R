#install.packages('RcppPlanc', repos = c('https://welch-lab.r-universe.dev', 'https://cloud.r-project.org'))
#devtools::install_github('linxihui/NNLM')
library(rliger)
library(NMF)
library(NNLM)
library(Matrix)
library(moltenNMF)
library(dplyr)
####
maxcor <- function(A, B){
  cor_0 <- cor(A,B)
  apply(cor_0, 1, max)
}

diagcor <- function(A,B){
  L = ncol(A)
  sapply(seq_len(L) , function(i)cor(A[,i], B[,i]))  
}


poisloss <- function(obs, fit){
  -mean(obs*log(fit) - fit -lgamma(obs + 1)) 
}

rmse <- function(obs, fit){
  sqrt(mean((as.matrix(obs) - fit)^2))
}

fit_mnmf <- function(out){
  Vhat_t <- moltenNMF:::meanV_array(out)
  Vhat_t[[1]]%*%t(Vhat_t[[2]])
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

zeroprob_norm <- function(mean, sd=1){
  pnorm(0, mean, sd)
}

basefun <- function(x){
  cbind(tanh(x),sin(x),cos(x),sign(x))  
}

softplus <- function(x){log1p(exp(x))}

relu <- function(x){
  x[x < 0] <- 0
  return(x)
}

row_norm <- function(V){
  sweep(V, 1, rowSums(V), "/")  
}

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
#####


set_simpar <- function(nrow, ncol, L, logmu=0){
  V1 <- matrix(rlnorm(L*nrow, meanlog = 0, sdlog = seq(1, 0.1, length.out=L)), nrow = nrow, ncol = L, byrow = TRUE)
  V2 <- matrix(rlnorm(L*ncol, logmu, 1), nrow = L, ncol = ncol)
  ord = order(apply(V1, 2, var), decreasing=TRUE)
  V1 = V1[,ord]
  V2 = V2[ord,]
  list(mu=V1%*%V2, V1=V1, V2=V2)
}


samplemat_norm <-function(tpar, obs_sigma){
  with(tpar,
       relu(matrix(rnorm(length(mu), mu, obs_sigma), nrow = nrow(mu), ncol = ncol(mu)))
  )
}

samplemat_pois <-function(tpar){
  with(tpar, 
       matrix(rpois(length(mu), mu), nrow = nrow(mu), ncol = ncol(mu))
  )
}

####

logmu = c(-5, -2, 0)
n_rows <- c(1000, 2000, 5000)
n_cols <- c(2000)
L_rank <- c(10L, 15L, 30L)
repl <- 10 
df <- expand.grid(n_rows=n_rows,n_cols=n_cols, logmu=logmu, L_rank=L_rank)
dim(df)

res_bm <- vector("list", nrow(df))
res_cor <- vector("list", nrow(df))
sparsity <- numeric(nrow(df))
pb <- txtProgressBar(0, nrow(df), style = 3)
for(i in 1:nrow(df)){
  set.seed(i)
  L <- df$L_rank[i]
  nrow = df$n_rows[i]
  ncol = df$n_cols[i]
  maxit=1000
tpar <- set_simpar(nrow, ncol, L, logmu = df$logmu[i])
#hist(tpar$mu, breaks = "FD")
#sp = with(tpar, mean(zeroprob_norm(mu, obs_sigma)))
sp = with(tpar, mean(exp(-tpar$mu)))
sparsity[i] <- sp
#Z = samplemat_norm(tpar, obs_sigma = 1)
Z = samplemat_pois(tpar)
Z = as(Z, "TsparseMatrix")

rank = L
t0 <- system.time2({
  decomp0 <- moltenNMF:::NMF2D_vb(Z, rank=rank, iter=maxit,
                                               display_progress = FALSE)
  })

t1 <- system.time2({
  decomp1 <- NMF::nmf(as.matrix(Z), rank = rank, maxIter=maxit, eps=0)
})

t2 = system.time2({
  decomp2 <- nnmf(as.matrix(Z), k = rank, rel.tol = 0, max.iter = maxit, loss = "mkl", verbose=0, method = "lee")
})

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
Vhat_t <- moltenNMF:::rearrange_cols(Vhat_t, FUN = var, normalize = FALSE, decreasing = TRUE)

fit_0 <- Vhat_t[[1]]%*%t(Vhat_t[[2]])
fit_1 <- NMF::basis(decomp1)%*%NMF::coef(decomp1)
fit_2 <- decomp2$W%*%decomp2$H
fit_3 <- fit_liger(liger_obj)


Vhat_1 = list(NMF::basis(decomp1), t(NMF::coef(decomp1)))
Vhat_1 <- moltenNMF:::rearrange_cols(Vhat_1, FUN = var, normalize = FALSE, decreasing = TRUE)

Vhat_2 <- list(decomp2$W, t(decomp2$H))
Vhat_2 <- moltenNMF:::rearrange_cols(Vhat_2, FUN = var, normalize = FALSE, decreasing = TRUE)

Vhat_3 <- list(extract_row_liger(liger_obj), t(getMatrix(liger_obj, slot = "H")[[1]]))
Vhat_3 <- moltenNMF:::rearrange_cols(Vhat_3, FUN = var, normalize = FALSE, decreasing = TRUE)


t0 = mutate(t0, poisloss = poisloss(as.matrix(Z), fit_0), 
             rmse = rmse(as.matrix(Z), fit_0), 
            rmse_t = rmse(tpar$mu, fit_0), 
            method="proposed")

t1 = mutate(t1, poisloss = poisloss(as.matrix(Z), fit_1),
             rmse = rmse(as.matrix(Z), fit_1), 
            rmse_t = rmse(tpar$mu, fit_1), 
            method="NMF")

t2 = mutate(t2, poisloss = poisloss(as.matrix(Z), fit_2),
             rmse = rmse(as.matrix(Z), fit_2),
            rmse_t = rmse(tpar$mu, fit_2), 
            method="NNLM")

t3 = mutate(t3, poisloss = poisloss(as.matrix(Z), fit_3),
            rmse = rmse(as.matrix(Z), fit_3),
            rmse_t = rmse(tpar$mu, fit_3), 
            method="rliger")

bm = rbind(t0,t1,t2,t3)

rho <- cbind(diagcor(Vhat_t[[1]], tpar$V1),diagcor(Vhat_1[[1]], tpar$V1), 
             diagcor(Vhat_2[[1]], tpar$V1),diagcor(Vhat_3[[1]], tpar$V1))

  res_bm[[i]] <- bm
  res_cor[[i]] <- rho
  setTxtProgressBar(pb,i)
}


save(res_bm, res_cor, sparsity, file = "vsNMF_batch.Rdata")


