sparse_cate <- function(x,repr="C"){
  if(!is.factor(x)){
    x <- factor(x)
    warning("auto-converted to factor")
  }
  xi <- as.integer(x)
  xi <- xi[!is.na(x)]
  xp <- seq(1,length(x))
  xp <- xp[!is.na(x)]
  val <- rep(TRUE, length(xi))
  m <- Matrix::sparseMatrix(i = xi, j = xp, x = val, repr=repr)
  rownames(m) <- levels(x)
  return(m)
}


sparse_onehot <- function(object,
                          data = environment(object),
                          #dummy_m = 0L,
                          xlev = NULL,
                          sep = "_",
                          interaction_operator = ":",
                          na.action='na.pass',
                          repr = "C"){
  data <- model.frame(object, data, xlev=xlev, na.action=na.action)
  t <- if(missing(data)) terms(object) else terms(object, data=data)
  labs <- attr(t, "term.labels")
  lx <- vector("list", length = length(labs)) 
  for(i in 1:length(labs)){
    if(any(grepl(interaction_operator, labs[i]))){
      lab2 <- unlist(strsplit(labs[i], interaction_operator))
      cmat <- sapply(lab2, function(x)as.character(data[[x]]))
      ct <- apply(cmat, 1, paste, collapse=":")
      lx[[i]] <- sparse_cate(ct, repr = repr)
    }else{
      lx[[i]] <- sparse_cate(data[[labs[i]]], repr = repr)      
    }
  }
  #if(dummy_m>=1L){
  #  m <- matrix(FALSE, dummy_m, nrow(data))
  #  m <- as(m, "lMatrix")
  #  rownames(m) <- 1:dummy_m
  #  lx[[length(labs)+1]] <- m
  #  labs <- c(labs, "dummy")
  #}
  X <- do.call("rbind", lx)
  X <- t(X)
  clabs <- unlist(sapply(lx, rownames))
  len <- unlist(sapply(lx, nrow))
  colnames(X) <- paste(rep(labs, len), clabs, sep=sep)
  attr(X, "indices") <- c(0L, cumsum(len))
  attr(X, "term.labels") <- labs
  attr(X, "value.labels") <- clabs
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
