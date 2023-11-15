matbarplot <- function(m, xgroup=NULL,
                       col=hcl.colors, ...){
  L <- ncol(m)
  cols <- col(L)
  bp <- barplot(t(m), col = cols,
                xaxt = "n", 
                ann = FALSE, axisnames = FALSE,
                border = cols, ...)
  uv <- unique(xgroup)
  nv <- length(uv)
  sta <- 1L
  en <- 0L
  if(!is.null(xgroup)){ # annotate x-axis
    for(i in 1:nv){
      en <- en + sum(xgroup==uv[i])
      axis(1, c(bp[sta], bp[en]), labels = rep("",2))
      axis(1, (bp[sta]+bp[en])*0.5, labels = uv[i], tick=FALSE)
      sta <- en+1L
    }
  }
}