library(moltenNMF)
library(Matrix)
library(dplyr)
library(bench)
library(ggplot2)
library(rliger)
library(readr)
# library(foreach)
# library(parallel)

poisloss <- function(obs, fit){
  -mean(obs*log(fit)-fit - lgamma(obs+1)) 
}

rmse <- function(obs, fit){
  sqrt(mean((as.matrix(obs) - fit)^2))
}

datapath <- scan("datapath.txt", what = character())

Y = readMM(datapath[4])
rownames(Y) <- paste0("gene", 1:nrow(Y))  # 遺伝子名
colnames(Y) <- paste0("cell", 1:ncol(Y))  # 細胞名
#dim(Y)
liger_obj = createLiger( list(sample1=Y) )
# liger_obj = rliger::normalize(liger_obj)
# liger_obj = rliger::selectGenes(liger_obj)
# liger_obj = rliger::scaleNotCenter(liger_obj)

liger_obj@datasets$sample1@normData <- liger_obj@datasets$sample1@rawData
liger_obj@datasets$sample1@scaleData <- liger_obj@datasets$sample1@normData
liger_obj@varFeatures <- rownames(liger_obj@datasets[[1]])
maxit = 10
bm2 = bench::mark({
  liger_obj <- runOnlineINMF(liger_obj,
                             k = 10L, nCores = 1,
                             minibatchSize = 1e+3L,
                             maxEpochs = maxit, 
                             verbose = FALSE)
}, iterations = 1)

print(bm2$total_time)

W = getMatrix(liger_obj, slot = "W")
V = getMatrix(liger_obj, slot = "V")
H = getMatrix(liger_obj, slot = "H")
fit_l = as.matrix((W + V$sample1)%*%H$sample1)
dim(fit_l)
bm2 = mutate(bm2, poisloss = poisloss(as.matrix(Y), fit_l),
             rmse = rmse(as.matrix(Y), fit_l))
#Y_l = getMatrix(liger_obj, "normData")

n_baches = 1e+3L
bm1 = bench::mark({
  out <- moltenNMF:::NMF2D_svb_t1(Y, rank = 10,
                             n_epochs = 10, n_baches = n_baches,
                             lr_param = c(15,0.9), lr_type = "exponential")
}, iterations = 1)

Vhat_t <- meanV_pois(out)
head.matrix(Vhat_t[[1]])
fit_0 <- Vhat_t[[1]]%*%t(Vhat_t[[2]])
bm1 = mutate(bm1, poisloss = poisloss(as.matrix(Y), fit_0),
            rmse = rmse(as.matrix(Y), fit_0))
bm1$total_time
bm2$total_time

# ggplot(data = NULL, aes(x=log1p(fit_l), y=as.matrix(log1p(Y_l$sample1))))+
#   geom_bin2d(aes(fill=after_stat(log1p(count))))+
#   theme_bw()