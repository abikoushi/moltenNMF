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
  out_b <- moltenNMF:::mNMF_bsvb(
    Y1,
    X = X1,
    N = nrow(X),
    L = L,
    n_epochs = 200,
    M_max = 10L,
    display_progress = TRUE
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
    M_max = 10,
    display_progress = TRUE
  )
})

#is.list(out_s$shape/out_s$rate)
plot(out_s$ELBO[-1], type = "l")

f = product_m(X, out_s$shape / out_s$rate)
plot(f, Y, pch = ".", col = rgb(0, 0, 0, 0.5))
abline(0, 1, col = "royalblue")

#TODO:batch shape update みなおし
