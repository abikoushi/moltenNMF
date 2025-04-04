library(moltenNMF)
library(Matrix)
library(ggplot2)
library(rbenchmark)

path <- scan("datapath.txt", what = character())
M = readRDS(path[1])
writeMM(M, file = path[2])
tpath <- path[2]
s = moltenNMF:::size_mtx(tpath)
system.time({
  resmv = moltenNMF:::rowmeanvar_mtx(n_row = s[1], n_col = s[2], readtxt = tpath, n_header = 2)
})
# ユーザ   システム       経過  
# 297.92       6.31     308.78 
saveRDS(resmv, file = "meanvar_MOCA.rds")

plot(resmv$mean, sqrt(resmv$var), col=rgb(0,0,0,0.1), cex=0.5, pch=16)
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
  moltenNMF:::rowfilter_mtx(readtxt = tpath, writetxt = "MOCA_2000.mtx", rowind = wch)  
})
# ユーザ   システム       経過  
# 1012.34      19.52    1043.49 

