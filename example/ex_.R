library(moltenNMF)
library(Matrix)
library(raster)
library(Rcpp)
library(jpeg)

datraw <- readJPEG("./example/gahag-002124.jpeg")
dat <- datraw[51:249,,] #trimming
plot(as.raster(dat))
ds <- dim(dat)
N <- length(dat)
#set.seed(1);ind <- sample.int(N, 0.5*N)
xall <- expand.grid(row=factor(1:ds[1]),
                    col=factor(1:ds[2]),
                    rgb=factor(1:ds[3]))
array(1:27,dim=c(3,3,3))
Y <- as.vector(dat)
tail(xall)

sX <- moltenNMF::sparse_onehot(~row+col+rgb, data=xall)

out <- mNMF_vb.default(1000*Y, sX, L=10, offset = rep(1000, length(Y)),  iter=200)
plot(out$ELBO[-1], type="l")

fit <- product_m(~row+col+rgb, data=xall, out$shape/out$rate)
plot(fit,Y,pch=".")
fit[fit>1] <- 1
fit[fit<0] <- 0
plot(as.raster(array(fit, dim=ds)))

