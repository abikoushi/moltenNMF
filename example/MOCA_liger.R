library(Matrix)
library(ggplot2)
library(rliger)
library(dplyr)

path <- scan("datapath.txt", what = character())
tpath <- path[2]
Y = readMM(tpath)
#dim(Y)
# log(0.95, base = 1-10^3/(2000*8337))
rownames(Y) <- paste0("gene", 1:nrow(Y))
colnames(Y) <- paste0("cell", 1:ncol(Y))
liger_obj = createLiger( list(sample1=Y) )
# liger_obj = rliger::normalize(liger_obj)
# liger_obj = rliger::selectGenes(liger_obj)
# liger_obj = rliger::scaleNotCenter(liger_obj)
liger_obj@datasets$sample1@normData <- liger_obj@datasets$sample1@rawData
liger_obj@datasets$sample1@scaleData <- liger_obj@datasets$sample1@normData
liger_obj@varFeatures <- rownames(liger_obj@datasets[[1]])
maxit = 100
rank = 30L
bm_l1 = bench::mark({
  liger_obj <- runOnlineINMF(liger_obj,
                             k = rank, nCores = 1,
                             minibatchSize = 1e+3L,
                             maxEpochs = maxit, 
                             verbose = FALSE)
}, iterations = 1)

save(bm_l1, liger_obj, file = "resMOCA_liger.Rdata")

maxit = 500
bm_l2 = bench::mark({
  liger_obj <- runOnlineINMF(liger_obj,
                             k = rank, nCores = 1,
                             minibatchSize = 1e+3L,
                             maxEpochs = maxit, 
                             verbose = FALSE)
}, iterations = 1)

save(bm_l2, liger_obj, file = "resMOCA_liger_2.Rdata")
