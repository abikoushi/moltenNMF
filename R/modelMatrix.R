sparse_cate <- function(x, fill){
  if(!is.factor(x)){
    x <- factor(x)
    warning("auto-converted to factor")
  }
  xi <- as.integer(x)
  xp <- seq(0L,length(x))
  m <- Matrix::sparseMatrix(i = xi, p = xp, x = rep(fill, length(x)))
  rownames(m) <- levels(x)
  return(m)
}

#' @export
sparse_onehot <- function(object,
                          data = environment(object),
                          xlev = NULL,
                          sep = "_",
                          na.action='na.pass',
                          fill = TRUE){
  data <- model.frame(object, data, xlev=xlev, na.action=na.action)
  t <- if(missing(data)) terms(object) else terms(object, data=data)
  labs <- attr(t, "term.labels")
  lx <- vector("list", length = length(labs))
  for(i in 1:length(labs)){
    lx[[i]] <- sparse_cate(data[[labs[i]]], fill = fill)
  }
  X <- do.call("rbind", lx)
  X <- t(X)
  clabs <- sapply(lx, rownames)
  len <- sapply(lx, nrow)
  colnames(X) <- paste(rep(labs, len), unlist(clabs), sep=sep)
  attr(X, "indices") <- c(0L, cumsum(len))
  attr(X, "term.labels") <- labs
  attr(X, "assign") <- rep(1:length(lx), len)
  return(X)
}

# sparse_model_matrix_b <- function(object,
#                                   data = environment(object),
#                                   xlev = NULL,
#                                   transpose = FALSE,
#                                   drop.unused.levels = FALSE,
#                                   row.names = TRUE,
#                                   sep = "_",
#                                   verbose = FALSE,
#                                   na.action='na.pass'){
#   data <- model.frame(object, data, xlev=xlev, na.action=na.action)
#   t <- if(missing(data)) terms(object) else terms(object, data=data)
#   X <- Matrix:::model.spmatrix(trms = t,
#                                mf = data,
#                                transpose = transpose,
#                                drop.unused.levels = drop.unused.levels,
#                                row.names = row.names,
#                                sep = sep,
#                                verbose = verbose)
#   m <- attr(X, "assign")
#   ind <- which(diff(m)!=0)
#   Y <- methods::as(X, "lgCMatrix") # equivalent to X > 0
#   attr(Y, "indices") <- c(0, ind, length(m))
#   attr(Y, "term.labels") <- attr(t, "term.labels")
#   dropK <- integer(length = length(ind)+1L)
#   dropK[-1] <- 1L
#   attr(Y, "dropK") <- dropK
#   attr(Y, "assign") <- attr(X, "assign")
#   return(Y)
# }
