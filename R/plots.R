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

matplot2 <- function(V,
                     xlab = "",
                     ylab = "",
                     palette = "black",
                     colour = NULL,
                     mar=c(5.1, 4.1, 4.1, 8.3),
                     legend = NULL, 
                     legendpos ="topright",
                     legendlwd = 2,
                     xshift = -0.24,
                     lwd=0.5, lty=1, 
                     type="l",...){
  if(is.null(colour)){
    col = palette
  }else{
    col = palette[colour]
  }
  if(is.null(legend) & !is.null(colour)){
    legend = levels(colour)
  }
  if(is.null(legend)){
    matplot(V, xaxt="n", xlab=xlab, ylab=ylab, type=type, lwd=lwd, col=col, lty=lty, ...)
  }else{
    oldpar = par(mar=mar, xpd=TRUE)
    matplot(V, xaxt="n", xlab=xlab, ylab=ylab, type=type, lwd=lwd, col=col, lty=lty, ...)
    legend(legendpos, inset=c(xshift,0), legend=legend, lty=lty, col=palette, lwd=legendlwd)
    par(oldpar)    
  }
  axis(side = 1, at = 1:nrow(V), labels = rownames(V))
}