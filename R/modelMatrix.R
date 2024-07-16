sparse_cate <- function(x, repr="C", binary_dummy=FALSE){
  if(!is.factor(x)){
    x <- factor(x)
    warning("auto-converted to factor")
  }
  if(binary_dummy & length(levels(x)) == 2L){
    xi <- as.integer(x)-1L
    f <- xi > 0L & !is.na(x)
    xi <- xi[f]
    xp <- seq(1,length(xi))
    m <- Matrix::sparseMatrix(i = xi, j = xp,
                              dims = c(1L, length(x)),
                              repr=repr)
    rownames(m) <- levels(x)[2]
  }else{
    xi <- as.integer(x)
    xi <- xi[!is.na(x)]
    xp <- seq(1,length(x))
    xp <- xp[!is.na(x)]
    m <- Matrix::sparseMatrix(i = xi, j = xp,
                              dims = c(length(levels(x)), length(x)),
                              repr=repr)
    rownames(m) <- levels(x)
  }
  return(m)
}


sparse_onehot <- function(object,
                          data = environment(object),
                          #dummy_m = 0L,
                          xlev = NULL,
                          sep = "_",
                          interaction_operator = ":",
                          na.action='na.pass',
                          repr = "C",
                          binary_dummy=FALSE){
  data <- model.frame(object, data, xlev=xlev, na.action=na.action)
  t <- if(missing(data)) terms(object) else terms(object, data=data)
  labs <- attr(t, "term.labels")
  lx <- vector("list", length = length(labs)) 
  for(i in 1:length(labs)){
    if(any(grepl(interaction_operator, labs[i]))){
      lab2 <- unlist(strsplit(labs[i], interaction_operator))
      cmat <- sapply(lab2, function(x){
        ifelse(is.na(data[[x]]),NA_character_,as.character(data[[x]]))
        })
      ct <- ifelse(apply(is.na(cmat),1,any), NA_character_, apply(cmat, 1, paste, collapse=":"))
      lx[[i]] <- sparse_cate(ct, repr = repr, binary_dummy = binary_dummy)
    }else{
      lx[[i]] <- sparse_cate(data[[labs[i]]], repr = repr, binary_dummy = binary_dummy)      
    }
  }
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
