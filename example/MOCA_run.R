library(moltenNMF)
library(Matrix)
library(ggplot2)
library(dplyr)
path <- scan("datapath.txt", what = character())
tpath <- path[2]
Y = readMM(tpath)
#log10(nnzero(Y))
#[1] 8.956141
maxit = 100
rank = 30L
set.seed(1234)
bm2 = bench::mark({
m_obj <- moltenNMF:::NMF2D_svb(Y, rank = rank,
                               n_epochs = maxit, n_baches = 1e+7L,
                               prior_shape = 1, prior_rate = 1,
                               lr_param = c(16, 0.9), lr_type = "exponential")
}, iterations = 1)

save(bm2, m_obj, file = "resMOCA.Rdata")
