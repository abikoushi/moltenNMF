library(ggplot2)
library(dplyr)
library(bench)
load("vsNMF_batch.Rdata")
logmu = c(-5, -2, 0)
n_rows <- c(1000, 2000, 5000)
n_cols <- c(2000)
L_rank <- c(10L, 15L, 30L)
repl <- 10 
df <- expand.grid(n_rows=n_rows,n_cols=n_cols, logmu=logmu, L_rank=L_rank)
dim(df)

df$L_rank
df$n_cols

t1 = sapply(res_bm, function(x)x$elapsed[1])
t2 = sapply(res_bm, function(x)x$elapsed[2])
t3 = sapply(res_bm, function(x)x$elapsed[3])
t4 = sapply(res_bm, function(x)x$elapsed[4])
plot(t1, type="l",col=1)
lines(t2, col=2)
lines(t3, col=3)
lines(t4, col=4)

t1 = sapply(res_bm, function(x)x$mem_used[1])
t2 = sapply(res_bm, function(x)x$mem_used[2])
t3 = sapply(res_bm, function(x)x$mem_used[3])
t4 = sapply(res_bm, function(x)x$mem_used[4])
plot(t1, type="l",col=1)
lines(t2, col=2)
lines(t3, col=3)
lines(t4, col=4)

s1 = sapply(res_bm, function(x)x$sparsity[1])
s2 = sapply(res_bm, function(x)x$sparsity[2])
s3 = sapply(res_bm, function(x)x$sparsity[3])
s4 = sapply(res_bm, function(x)x$sparsity[4])

plot(sparsity)

matplot(res_cor[[1]], type = "l")

rhom = t(sapply(res_cor, colMeans))
matplot(rhom, type = "l")
