library(moltenNMF)
library(Matrix)
library(tidyr)
library(dplyr)
library(ggplot2)

set_attr_modelmat <- function(X) {
  if (!is.null(attr(X, "assign"))) {
    attr(X, "indices") <- c(0L, cumsum(rle(attr(X, "assign"))$lengths))
  }
  labs <- names(attr(X, "contrasts"))
  if (!is.null(labs)) {
    attr(X, "term.labels") <- labs
    if (!is.null(X@Dimnames[[2]])) {
      attr(X, "value.labels") <- sub(
        paste(labs, collapse = "|"),
        "",
        X@Dimnames[[2]]
      )
    }
  }
  return(X)
}

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
V <- matrix(rlnorm(L * D, -1, 1), D, L)
ord = order(apply(log(V), 2, var), decreasing = TRUE)
V = V[, ord]
lambda <- product_m.default(X, V)
Y <- rpois(N, lambda)

st = system.time({
  out_d <- moltenNMF::mNMF_vb.default(
    Y,
    X = X,
    L = L,
    iter = 500,
    display_progress = TRUE
  )
})
st

Xs <- sparse.model.matrix(~., data = df)
Xs = set_attr_modelmat(Xs)

wch = which(Y > 0)
Y_sp = sparseVector(Y[wch], wch, length = length(Y))
st2 = system.time({
  out_sb <- moltenNMF::mNMF_vb.default(
    Y_sp,
    X = Xs,
    L = L,
    iter = 500,
    display_progress = TRUE
  )
})

st2 / st


Xs1 = slice_rows(Xs, wch)
st_bs <- system.time({
  out_bs <- moltenNMF:::mNMF_bsvb(
    Y_sp@x,
    X = Xs1,
    N = nrow(X),
    L = L,
    n_epochs = 500,
    M_max = 10L,
    display_progress = TRUE
  )
})

st_bs

plot(out_d$ELBO[-1], type = "l")
plot(out_bs$ELBO[-1], type = "l", col = "royalblue", lty = 2)
plot(out_bs$ELBO[-1], type = "l")

V_bs <- out_bs$shape / out_bs$rate
f_bs = product_m(Xs, V_bs)
V_d <- out_d$shape / out_d$rate
f_d <- moltenNMF::product_m(X, V_d)
V_sb <- out_sb$shape / out_sb$rate
f_sb <- moltenNMF::product_m(Xs, V_sb)

ggplot() +
  geom_point(aes(x = f_d, y = as.matrix(Y)), alpha = 0.25, shape = 1) +
  geom_point(
    aes(x = f_bs, y = as.matrix(Y)),
    alpha = 0.25,
    colour = "royalblue",
    shape = 2
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "lightgrey") +
  theme_bw()


ba05 = rgb(0, 0, 0, 0.5)
reV_s = rearrange_winner_ord(V, V_sb)
plot(log(V), log(reV_s$V), pch = 1, col = ba05)
abline(0, 1, col = "lightgrey")

reV_d = rearrange_winner_ord(V, V_d)
plot(log(V), log(reV_d$V), col = ba05)
abline(0, 1, col = "lightgrey")
