library(moltenNMF)
library(Matrix)
library(dplyr)
library(bench)
# library(foreach)
# library(parallel)
library(ggplot2)
library(rliger)
library(readr)

# X = sparseMatrix(i=1:4, j=c(1,2,2,4), x=1:4)
# X@p
# #extract a column
# X@x[X@i[(X@p[2]+1):X@p[3]]+1]

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
                             k=10L, nCores = 1,
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
Y_l = getMatrix(liger_obj, "normData")
ggplot(data = NULL, aes(x=log1p(fit_l), y=as.matrix(log1p(Y_l$sample1))))+
  geom_bin2d(aes(fill=after_stat(log1p(count))))+
  theme_bw()

#dim(Y)
## [1]  33514 120961
#nnzero(Y)
##[1] 239634370
#print(nnzero(Y)/length(Y))
## [1] 0.05911225
wch = which(rowSums(Y) > 0)
length(wch)

#log10(nnzero(Y))
##[1] 8.379549

out <- moltenNMF:::NMF2D_svb(Y, rank = 10,
                             n_epochs = 50, n_baches = 1e+7L,
                             lr_param = c(1.5,0.6), lr_type = "power")
