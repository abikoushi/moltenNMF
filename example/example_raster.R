#install.packages("devtools")
library(raster)
library(jpeg)
library(moltenNMF)
library(animation)

mat4img <- function(x){
  t(apply(x, 2, rev))
}

dat <- readJPEG("./example/Shunkosai_Hokuei_Obake.jpeg")
plot(as.raster(dat))
ds <- dim(dat)
N <- prod(ds)
set.seed(1);ind <- sample.int(N, 0.9*N)
dat2 <- dat
dat2[ind] <- 1
plot(as.raster(dat2))
x <- expand.grid(row=factor(1:ds[1]),
                 col=factor(1:ds[2]),
                 rgb=factor(1:ds[3]))
coef <- 100
Y <- floor(as.vector(coef*dat))[-ind]
X <- x[-ind,]
system.time({
  out <- mNMF_vb(Y~row+col+rgb-1, data=X, L=30, iter=500)  
})
# user  system elapsed 
# 83.577   7.809 101.733 
plot(out$ELBO, type="l")

# dx <- sparse_model_matrix_b(~row+col+rgba-1,data = x)
# y_pred <- rpredictor_mNMF(dx, 100, out$shape, out$rate)

fit <- product_m(~ row+col+rgb-1, data=x, out$shape/out$rate)
Y2 <- as.vector(dat)
Y2[ind] <- 0

f2 <- fit/coef
plot(as.raster(array(Y2, dim=ds)))
f2[f2>1] <- 1
Y2imp <- Y2
Y2imp[ind] <- f2[ind]
plot(as.raster(array(Y2imp, dim=ds)))

jpeg("hokusai_oiwa.jpeg",width = 480, height = 240)
par(mfrow=c(1,3), mai=c(0,0,0,0))
plot(as.raster(dat))
plot(as.raster(array(Y2, dim=ds)))
plot(as.raster(array(Y2imp, dim=ds)))
dev.off()

