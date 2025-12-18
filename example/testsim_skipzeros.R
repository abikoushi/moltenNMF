library(moltenNMF)
library(Matrix)
library(tidyr)
library(dplyr)
library(ggplot2)

rearrange_winner_ord <- function(V, V_s) {
  cmat = cor(V, V_s)
  ord = integer(ncol(V))
  ord[1] = which.max(cmat[1, ])
  for (i in 2:L) {
    ord[i] = which(cmat[i, ] == max(cmat[i, -ord[1:(i - 1)]]))
  }
  list(V = V_s[, ord], cor = diag(cmat[, ord]))
}

i = 1
L <- 5L
ncols = c(100, 500, 1000)
df <- as.data.frame(expand.grid(
  row = factor(1:100),
  col = factor(1:ncols[i]),
  depth = factor(1:2)
))

X <- sparse_onehot(~., data = df)
N <- nrow(X)
D <- ncol(X)
set.seed(123456789)
V <- matrix(rlnorm(L * D, 0, 1), D, L)
ord = order(apply(log(V), 2, var), decreasing = TRUE)
V = V[, ord]
lambda <- product_m.default(X, V)
Y <- rpois(N, lambda)
wch = which(Y > 0)

Y1 = Y[wch]
X1 = slice_rows(X, wch)
length(wch) / length(Y)

system.time({
  out_b <- moltenNMF:::mNMF_svb_batch(
    Y1,
    X = X1,
    N = nrow(X),
    L = L,
    n_epochs = 100,
    lr_param = c(1, 0.8),
    lr_type = "exponential",
    M_max = 10L,
    display_progress = FALSE
  )
})
plot(out_b$ELBO[-1], type = "l")

f = product_m(X, out_b$shape / out_b$rate)
plot(f, Y, pch = ".", col = rgb(0, 0, 0, 0.5))
abline(0, 1, col = "royalblue")

system.time({
  out_s <- moltenNMF:::mNMF_svb(
    Y1,
    X = X1,
    N = nrow(X),
    L = L,
    n_batches = 1000,
    n_epochs = 200,
    lr_param = c(1, 0.8),
    lr_type = "exponential",
    M_max = 0,
    display_progress = TRUE
  )
})

#is.list(out_s$shape/out_s$rate)
plot(out_s$ELBO[-1], type = "l")

f = product_m(X, out_s$shape / out_s$rate)
plot(f, Y, pch = ".", col = rgb(0, 0, 0, 0.5))
abline(0, 1, col = "royalblue")


####
#batch
####

system.time({
  out_d <- moltenNMF::mNMF_vb.default(
    Y,
    X = X,
    L = L,
    iter = 500,
    a = 1,
    b = 1,
    display_progress = TRUE
  )
})


wch = which(Y > 0)
Y_sp = sparseVector(Y1, wch, length = length(Y))
system.time({
  out_sb = moltenNMF:::mNMF_vb.default(Y_sp, X, L = L, iter = 100)
})

plot(out_d$ELBO[-1], type = "l")
lines(out_sb$ELBO[-1], type = "l", col = "royalblue", lty = 2)


V_d <- out_d$shape / out_d$rate
f_d <- moltenNMF::product_m(X, V_d)
V_sb <- out_sb$shape / out_sb$rate
f_sb <- moltenNMF::product_m(X, V_sb)

ggplot() +
  geom_point(aes(x = f_d, y = as.matrix(Y)), alpha = 0.25, shape = 1) +
  geom_point(
    aes(x = f_sb, y = as.matrix(Y)),
    alpha = 0.25,
    colour = "royalblue",
    shape = 2
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "lightgrey") +
  theme_bw()

rearrange_winner_ord <- function(V, V_s) {
  cmat = cor(V, V_s)
  ord = integer(ncol(V))
  ord[1] = which.max(cmat[1, ])
  for (i in 2:L) {
    ord[i] = which(cmat[i, ] == max(cmat[i, -ord[1:(i - 1)]]))
  }
  list(V = V_s[, ord], cor = diag(cmat[, ord]))
}

ba05 = rgb(0, 0, 0, 0.5)
reV_s = rearrange_winner_ord(V, V_s)
plot(log(V), log(reV_s$V), pch = 1, col = ba05)
abline(0, 1, col = "lightgrey")

reV_d = rearrange_winner_ord(V, V_d)
plot(log(V), log(reV_d$V), col = ba05)
abline(0, 1, col = "lightgrey")

rownames(V_s)
