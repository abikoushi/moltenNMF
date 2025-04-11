\name{matplot2}
\alias{matplot2}
\title{Customized Matrix Plot with Optional Legend}
\description{
  A wrapper around \code{matplot} to visualize matrices with customizable colors, line styles, and an optional legend.
}
\usage{
matplot2(V,
         xlab = "",
         ylab = "",
         palette = "black",
         colour = NULL,
         mar = c(5.1, 4.1, 4.1, 8.3),
         legend = NULL,
         legendpos = "topright",
         legendlwd = 2,
         xshift = -0.24,
         lwd = 0.5,
         lty = 1,
         type = "l",
         ...)
}
\arguments{
  \item{V}{A numeric matrix to plot. Each column is plotted as a separate line.}
  \item{xlab}{Label for the x-axis.}
  \item{ylab}{Label for the y-axis.}
  \item{palette}{A color palette to use. Either a single color or a vector of colors.}
  \item{colour}{A factor or vector used to select colors from the palette.}
  \item{mar}{A numeric vector specifying the margins, passed to \code{par(mar = ...)}.}
  \item{legend}{Character vector of legend labels. If \code{NULL}, and \code{colour} is provided, levels of \code{colour} are used.}
  \item{legendpos}{Position of the legend (e.g., "topright", "bottomleft").}
  \item{legendlwd}{Line width used for legend symbols.}
  \item{xshift}{Horizontal shift for the legend inset.}
  \item{lwd}{Line width for the plot.}
  \item{lty}{Line type.}
  \item{type}{Plot type (e.g., "l" for lines, "p" for points).}
  \item{...}{Additional graphical parameters passed to \code{matplot}.}
}
\details{
  The function enhances the default \code{matplot} behavior by allowing grouped color schemes and optional legends. Useful for plotting factorized matrices (e.g., NMF components).
}
\value{
  No return value. Called for its side effect (a plot).
}
\seealso{
  \code{\link{matplot}}, \code{\link{legend}}, \code{\link{par}}
}
\examples{
set.seed(1)
V <- matrix(runif(50), nrow=10)
rownames(V) <- paste0("T", 1:10)
matplot2(V, palette=rainbow(5), colour=rep(1:5, each=2))
}
