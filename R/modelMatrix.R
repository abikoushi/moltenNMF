sparse_cate <- function(x, repr = "C", binary_dummy=FALSE){
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
  
  li <- vector("list", length = length(labs))
  lp <- vector("list", length = length(labs))
  l_labs <- vector("list", length = length(labs))
  len = integer(length(labs))
  last_p <- 0L
  for(i in 1:length(labs)){
    if(any(grepl(interaction_operator, labs[i]))){
      lab2 <- unlist(strsplit(labs[i], interaction_operator))
      cmat <- sapply(lab2, function(x){
        ifelse(is.na(data[[x]]), NA_character_, as.character(data[[x]]))
      })
      ct <- ifelse(apply(is.na(cmat),1,any), NA_character_, apply(cmat, 1, paste, collapse=":"))
      X_i <- t(moltenNMF:::sparse_cate(ct, repr = repr, binary_dummy = binary_dummy))
    }else{
      X_i <- t(moltenNMF:::sparse_cate(data[[labs[i]]], repr = repr, binary_dummy = binary_dummy))
    }
    l_labs[[i]] <- colnames(X_i)
    len[i] <- ncol(X_i)
    li[[i]] <- X_i@i
    lp[[i]] <- X_i@p[-1] + last_p
    last_p <- lp[[i]][length(lp[[i]])]
  }
  X = sparseMatrix(i = unlist(li) + 1L,  p = c(0L, unlist(lp)))
  clabs <- unlist(l_labs)
  # len <- unlist(sapply(lx, nrow))
  colnames(X) <- paste(rep(labs, len), clabs, sep=sep)
  attr(X, "indices") <- c(0L, cumsum(len))
  attr(X, "term.labels") <- labs
  attr(X, "value.labels") <- clabs
  attr(X, "assign") <- rep(1:length(labs), len)
  return( X )
}

# sparse_onehot <- function(object,
#                           data = environment(object),
#                           #dummy_m = 0L,
#                           xlev = NULL,
#                           sep = "_",
#                           interaction_operator = ":",
#                           na.action='na.pass',
#                           repr = "C",
#                           binary_dummy=FALSE){
#   data <- model.frame(object, data, xlev=xlev, na.action=na.action)
#   t <- if(missing(data)) terms(object) else terms(object, data=data)
#   labs <- attr(t, "term.labels")
#   lx <- vector("list", length = length(labs)) 
#   for(i in 1:length(labs)){
#     if(any(grepl(interaction_operator, labs[i]))){
#       lab2 <- unlist(strsplit(labs[i], interaction_operator))
#       cmat <- sapply(lab2, function(x){
#         ifelse(is.na(data[[x]]),NA_character_,as.character(data[[x]]))
#         })
#       ct <- ifelse(apply(is.na(cmat),1,any), NA_character_, apply(cmat, 1, paste, collapse=":"))
#       lx[[i]] <- sparse_cate(ct, repr = repr, binary_dummy = binary_dummy)
#     }else{
#       lx[[i]] <- sparse_cate(data[[labs[i]]], repr = repr, binary_dummy = binary_dummy)      
#     }
#   }
#   X <- do.call("rbind", lx)
#   X <- t(X)
#   clabs <- unlist(sapply(lx, rownames))
#   len <- unlist(sapply(lx, nrow))
#   colnames(X) <- paste(rep(labs, len), clabs, sep=sep)
#   attr(X, "indices") <- c(0L, cumsum(len))
#   attr(X, "term.labels") <- labs
#   attr(X, "value.labels") <- clabs
#   attr(X, "assign") <- rep(1:length(lx), len)
#   return(X)
# }

slice_rows <- function(x, i){
  out = x[i,]
  if(!is.null(attr(x,"indices"))){
    attr(out, "indices") = attr(x,"indices")    
  }
  if(!is.null(attr(x,"term.labels"))){
    attr(out, "term.labels") = attr(x, "term.labels")
  }
  if(!is.null(attr(x, "value.labels"))){
    attr(out, "value.labels") = attr(x, "value.labels")[i]
  }
  if(!is.null(attr(x, "assign"))){
    attr(out, "assign") = attr(x, "assign")
  }
  return(out)
}

append_new <- function(x, new, onehot = NULL, newterm = NULL){
  out = cbind(x, new)
  y_len = ncol(new)
  if(is.null(onehot)){
    onehot = all(rowSums(new) <= 1L)    
  }
  if(is.null(newterm)){
    newterm = paste(as.character(substitute(new)), collapse = "")
    #print(newterm)
  }
  ###
  ind = attr(x, "indices")
  if(!is.null(ind)){
    if(onehot){
      attr(out, "indices") = c(ind, ind[length(ind)] + y_len)
    }else{
      attr(out, "indices") = c(ind, ind[length(ind)] + 1:y_len)
    }
  }
  assig = attr(x, "assign")
  if(!is.null(assig)){
    if(onehot){
      attr(out, "assign") <- c(assig, rep(assig[length(assig)]+1L, y_len))
    }else{
      attr(out, "assign") <- c(assig, assig[length(assig)]+1L:y_len)
    }
  }
  if(!is.null(attr(x, "term.labels"))){
    attr(out, "term.labels") = c(attr(x, "term.labels"), newterm)
  }
  if(!is.null(attr(x, "value.labels"))){
    attr(out, "value.labels") = c(attr(x, "value.labels"), rownames(new))
  }
  colnames(out) <- c(colnames(x), paste0(newterm, 1L:y_len))
  return(out)
}
