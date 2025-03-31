library(moltenNMF)
library(Matrix)
library(dplyr)
library(bench)
library(foreach)
library(parallel)
library(rliger)
library(readr)
library(ggplot2)
#library(snow)
#library(Rmpi)
#cl <- makeCluster(10, type='MPI')
#registerDoParallel(cl)
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

set.seed(5147)

#bmres = foreach(i = 1:27, 
#        .packages = c("moltenNMF2", "Matrix","dplyr"," pryr")) %do% {
poisloss <- function(obs, fit){
  -mean(obs*log(fit)-fit - lgamma(obs+1)) 
}
rmse <- function(obs, fit){
  sqrt(mean((as.matrix(obs) - fit)^2))
}

condf = expand.grid(rank = c(10, 20, 30),
                    n_b = c(1e4L, 1e5L, 1e6L), 
                    forgetting = c(0.7, 0.8, 0.9))
path <- scan("datapath.txt", what = character())
Y = Matrix::readMM(path[11])
#dim(Y)
maxit = 10

Cells = read_csv(path[9])
Genes = read_csv(path[10], col_names = "Genes")

genepos = scan(path[12])

rownames(Y) <- Genes$Genes[genepos]
colnames(Y) <- Cells$cell_name
liger_obj = createLiger( list(sample1=as.matrix(Y)) )
# liger_obj = rliger::normalize(liger_obj)
# liger_obj = rliger::selectGenes(liger_obj)
# liger_obj = rliger::scaleNotCenter(liger_obj)

liger_obj@datasets$sample1@normData <- liger_obj@datasets$sample1@rawData
liger_obj@datasets$sample1@scaleData <- liger_obj@datasets$sample1@normData
liger_obj@varFeatures <- rownames(liger_obj@datasets[[1]])
condf$rank
i=1
bm2 = bench::mark({
  liger_obj <- runOnlineINMF(liger_obj,
                             k=condf$rank[i], nCores = 1,
                             minibatchSize = 1e+3L,
                             maxEpochs = maxit, 
                             verbose = FALSE)
}, iterations = 1)

print(bm2$total_time)

W = getMatrix(liger_obj, slot = "W")
V = getMatrix(liger_obj, slot = "V")
H = getMatrix(liger_obj, slot = "H")

tied_min<- function(value){
  mv = min(value)
  list(min = mv, ratio=mean(value == mv))
}

print( tied_min(W) )
print( tied_min(V$sample1) )
print( tied_min(H$sample1) )

fit_l = as.matrix((W + V$sample1)%*%H$sample1)
bm2 = mutate(bm2, poisloss = poisloss(as.matrix(Y), fit_l),
            rmse = rmse(as.matrix(Y), fit_l))
bm2

V1 = matrix(runif(nrow(Y)*10, 0, 0.1), nrow(Y), 10)
V2 = t(unname(apply(Y, 2, sample, size=10)))+0.0001
Vini = list(V1,V2)
n_baches = 1e+3L
bm = bench::mark({
  out_t <- moltenNMF:::NMF2D_svb_t1(Y, rank=condf$rank[i],
                        n_epochs = 10,
                        n_baches = n_baches,
                        Vini = Vini,
                        # lr_param = c(16, 0.9),
                        # lr_type = "exponential",
                        lr_param = c(0.1),
                        lr_type = "const",
                        display_progress = FALSE)
}, iterations = 1)

plot(out_t$ELBO, type = "l")
Vhat_t <- meanV_pois(out_t)
head.matrix(Vhat_t[[1]])
fit_0 <- Vhat_t[[1]]%*%t(Vhat_t[[2]])
bm = mutate(bm, poisloss = poisloss(as.matrix(Y), fit_0),
            rmse = rmse(as.matrix(Y), fit_0))
bm$poisloss
bm2$poisloss

bw = c(0.25, 0.25)
ggplot(data=NULL, aes(x=log1p(c(fit_l)), y=c(as.matrix(log1p(Y)))))+
  geom_bin2d(aes(fill=after_stat(log10(count))), binwidth = bw)+
  geom_abline(slope=1, intercept = 0, linetype=2)+theme_bw()

ggplot(data=NULL, aes(x=log1p(c(fit_0)), y=c(as.matrix(log1p(Y)))))+
  geom_bin2d(aes(fill=after_stat(log10(count))), binwidth = bw)+
  geom_abline(slope=1, intercept = 0, linetype=2)+theme_bw()

# bw = c(10, 10)
# ggplot(data=NULL, aes(x=fit_l, y=as.matrix(Y)))+
#   geom_bin2d(aes(fill=after_stat(log10(count))), binwidth = bw)+
#   geom_abline(slope=1, intercept = 0, linetype=2)+theme_bw()
# 
# ggplot(data=NULL, aes(x=fit_0, y=as.matrix(Y)))+
#   geom_bin2d(aes(fill=after_stat(log10(count))), binwidth=bw)+
#   geom_abline(slope=1, intercept = 0, linetype=2)+
#   theme_bw()
# 

#saveRDS(bmres, file = paste0("bm_Song_svb_",rank,".rds"))

