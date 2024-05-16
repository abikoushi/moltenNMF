grepV <- function(V, x,
                  normalize = FALSE,
                  simplifynames = FALSE){
  out <- V[grep(x,rownames(V)),]
  if(normalize){
    out <- sweep(out,2,rowSums(out),"/")
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
