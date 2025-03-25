library(moltenNMF)
library(Matrix)
library(dplyr)
library(ggplot2)

set_data_mf <- function(L, nrow, ncol, mu=0){
  W <- matrix(rnorm(nrow*L,0,1), ncol=L)
  H <- matrix(rnorm(L*ncol,0,1), nrow=L)
  W <- sweep(W,1,rowMeans(W)-mu)
  H <- sweep(H,2,colMeans(H)-mu)
  Y <- matrix(rpois(nrow*ncol, exp(W)%*%exp(H)), nrow, ncol)
  Y <- as(Y, "TsparseMatrix")
  list(Y=Y, trueW=W, trueH=H)
}

dat <- set_data_mf(3, 111, 133)
hist(as.matrix(dat$Y))
out <- moltenNMF:::NMF2D_vb(dat$Y, rank = 3, iter = 1000)
plot(out$ELBO[-1], type = "l")
V = moltenNMF:::meanV_array(out)
head(out$shape[[1]])
head(out$shape[[2]])
out$rate

ggplot(data = NULL, aes(x=c(V[[1]]%*%t(V[[2]])), y=c(as.matrix(dat$Y))))+
  geom_bin2d(aes(fill = after_stat(log10(count))))+
  geom_abline(intercept = 0, slope = 1, linetype=2, colour="grey")

nnzero(dat$Y)
#lr = 1000/nnzero(dat$Y)
out2 <- moltenNMF:::NMF2D_svb(dat$Y, rank = 3,
                             n_epochs = 100, n_baches = as.integer(2000),
                             prior_shape = 1, prior_rate = 1,
                             lr_param = c(16,0.8), lr_type = "exponential")
out2$ELBO
plot(out2$ELBO, type="l")
head(out$shape[[1]])
head(out$shape[[2]])
head(out2$shape[[1]])
head(out2$shape[[2]])
out$rate
out2$rate

V = moltenNMF:::meanV_array(out2)
V = moltenNMF:::rearrange_cols(V, normalize = FALSE)
moltenNMF:::matplot2(V = t(V[[1]]))

fit2 = V[[1]]%*%t(V[[2]])
ax_lim = c(0, max(max(dat$Y), max(fit2)))

ggplot(data = NULL, aes(x=c(fit2), y=c(as.matrix(dat$Y))))+
  geom_bin2d(aes(fill = after_stat(log10(count))))+
  geom_abline(intercept = 0, slope = 1, linetype=2, colour="grey")

