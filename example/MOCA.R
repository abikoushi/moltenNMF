library(moltenNMF)
library(Matrix)
library(ggplot2)
library(rbenchmark)
path <- scan("datapath.txt", what = character())

# tpath <- path[2]
# s = size_mtx(tpath)
# system.time({
#   resmv = rowmeanvar_mtx(n_row = s[1], n_col = s[2], readtxt = tpath, n_header = 2)
# })
#    user   system  elapsed 
#1315.509    7.882 1326.068 
#saveRDS(resmv, file = "meanvar_MOCA.rds")

#plot(resmv$mean, sqrt(resmv$var), col=rgb(0,0,0,0.1), cex=0.5, pch=16)
resmv <- readRDS("meanvar_MOCA.rds")

wch = which(resmv$var>0)
# which(!resmv$var>0)
# [1] 22571
#length(wch)
#[1] 26182
#length(resmv$var)
#26183


wch = order(resmv$var ,decreasing = TRUE)[1:2000]
cat(wch, file = "rowposMOCA.txt")

system.time({
  rowfilter_mtx(readtxt = tpath, writetxt = "MOCA_2000.mtx", rowind = wch)  
})
#    user   system  elapsed 
#3687.166   21.206 3716.336 

