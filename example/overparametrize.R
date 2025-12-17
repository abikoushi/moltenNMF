library(moltenNMF)
library(Matrix)

set_attr_modelmat <- function(X){
  if(!is.null(attr(X, "assign"))){
    attr(X, "indices") <- c(0L, cumsum(rle(attr(X, "assign"))$lengths))    
  }
  labs <- names(attr(X, "contrasts"))
  if(!is.null(labs)){
    attr(X, "term.labels") <- labs
    if(!is.null(X@Dimnames[[2]])){
      attr(X, "value.labels") <- sub(paste(labs, collapse = "|"), "", X@Dimnames[[2]])    
    }
  }
  return(X)
}

df = as.data.frame(Titanic)
f = Freq ~ Class + Sex + Age + Survived
X_s = sparse_onehot(f, data = df)

L=2
system.time({
  out_s <- mNMF_vb.default(df$Freq, X = X_s, L = L, iter = 1000)
})

str(X)

X2 = sparse.model.matrix(f, data = df)
X2 = set_attr_modelmat(X2)
str(X2)

system.time({
  out <- mNMF_vb.default(df$Freq, X = X2, L = L, iter = 1000)
})

plot(out_s$ELBO[-1], type="l", lty=2)
lines(out$ELBO[-1], type="l")

fit = product_m.default(X2, out$shape/out$rate)
fit_s = product_m.default(X_s, out_s$shape/out_s$rate)
plot(df$Freq, fit)
points(df$Freq, fit_s, col="royalblue", pch=2)
