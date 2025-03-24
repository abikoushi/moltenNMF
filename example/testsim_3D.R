library(rbenchmark)
library(moltenNMF)
library(Matrix)
library(dplyr)
library(ggplot2)

ncol=6
df = expand.grid(row=1:100, col=1:ncol, depth=1:3)
dim(df)
trueV <- list(matrix(rexp(100*2, 0.5), 100, 2),
              matrix(rexp(ncol*2, 0.5), ncol, 2),
              matrix(rexp(3*2, 0.5), 3, 2))
truef = moltenNMF:::product_array(trueV, as.matrix(df))
y = rpois(nrow(df), truef)

plot(truef, y)
abline(0, 1, col="lightgrey")

iter = 1000

y0 = y[y>0]
X0 = slice_rows(as.matrix(df), y>0)
system.time({
  out1 = moltenNMF:::NTF_vb(Y = y0,
                            X = X0,
                            rank = 2,
                            iter = iter,
                            dims = c(100, ncol, 3),
                            #weight = list(rep(1,100), rep(1,10), rep(1,3)),
                            prior_shape = 1, prior_rate = 1,
                            display_progress = TRUE)
})
# user  system elapsed 
# 0.365   0.007   0.378 

plot(out1$logprob, type = "l")
Vhat = moltenNMF:::meanV_array(out1)
fhat1 = moltenNMF:::product_array(V = Vhat, X = as.matrix(df))

ggplot(df, aes(x=fhat1, y=y))+
  geom_bin2d(aes(fill = after_stat(log10(count))))+
  geom_abline(intercept = 0, slope = 1, linetype=2)

###
length(y0)
head(X0)
system.time({
  out2 = moltenNMF:::NTF_svb(Y = y0,
                             X = X0,
                             rank = 2,
                             n_epochs = 20,
                             n_baches = 1000,
                             lr_param = c(16, 0.8),
                             dims = c(100, ncol, 3),
                             prior_shape = 1, prior_rate = 1,
                             display_progress = TRUE)
})

out2$rate
plot(out2$logprob, type = "l")
Vhat = moltenNMF:::meanV_array(out2)
fhat2 = moltenNMF:::product_array(V = Vhat, X = as.matrix(df))

ggplot(df, aes(x=fhat2, y=y))+
  geom_bin2d(aes(fill = after_stat(log10(count))))+
  geom_abline(intercept = 0, slope = 1, linetype=2)

