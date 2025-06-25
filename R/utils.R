grepV <- function(V, x,
                  normalize = FALSE,
                  simplifynames = FALSE){
  out <- V[grep(x,rownames(V)),,drop=FALSE]
  if(normalize){
    out <- sweep(out, 1, rowSums(out),"/")
  }
  if(simplifynames){
    row.names(out) <- gsub(x,"",row.names(out))
  }
  return(out)
}

simmilar_top_n <- function(V, ref,
                           n=10L,
                           decreasing=FALSE){
  d <- rowSums(sweep(V, 2, ref)^2)
  ord <- order(d, decreasing = decreasing)[1:n]
  out <- data.frame(V[ord,], dist=d[ord])
  return(out)
}

learning_rate <- function(x, lr_param, lr_type){
  sapply(x, check_lr, lr_param = lr_param, lr_type = lr_type)  
}

size_mtx <- function(file_path){
  con <- file(file_path, open = "r") #Open for reading in text mode
  size <- scan(con, what=integer(), comment.char = "%", nlines = 1, skip = 1, quiet = TRUE)
  close(con)
  size
}

rowmeanvar_mtx <- function(file_path, n_header=2){
  size = size_mtx(file_path = file_path)
  rowmeanvar_txt(n_row = size[1],
                 n_col = size[2],
                 file_path,
                 n_header = n_header)
}

writebin_spmat <- function(object, filepath_x, filepath_y,
                           index_decrement=0L){
  if(all(class(object) != "dgTMatrix")){
    object = as(object, "TsparseMatrix")    
  }
  writeBin(object@x, filepath_y)
  writeBin(c(rbind(object@i-index_decrement, object@j-index_decrement)), filepath_x)
}
