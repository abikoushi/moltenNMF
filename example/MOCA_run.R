library(moltenNMF)
library(Matrix)
library(ggplot2)
library(rbenchmark)
library(rliger)
library(dplyr)
poisloss <- function(obs, fit){
  -mean(obs*log(fit)-fit - lgamma(obs+1)) 
}

rmse <- function(obs, fit){
  sqrt(mean((as.matrix(obs) - fit)^2))
}

tpath <- path[13]
Y = readMM(tpath)

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


W = getMatrix(liger_obj, slot = "W")
V = getMatrix(liger_obj, slot = "V")
H = getMatrix(liger_obj, slot = "H")
fit_l = as.matrix((W + V$sample1)%*%H$sample1)

writeBin(c(fit_l), con = "fitliger.bin")

log(0.95, base = 1-10^3/(2000*8337))
