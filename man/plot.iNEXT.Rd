\name{plot.iNEXT}
\alias{plot.iNEXT}
\title{Plot Method for an iNEXT Object}
\usage{
  plot.iNEXT(x, style = "N2D", col = 1, ...)
}
\arguments{
  \item{x}{a \code{\link{iNEXT}} objext computed by
  \code{\link{iNEXT.Ind}} or \code{\link{iNEXT.Sam}}}

  \item{style}{three different plotting style = c("N2D",
  "N2SC", "SC2D"),}

  \item{col}{a specification for the default plotting
  color, see \code{\link{par}} for detail.}

  \item{...}{further plotting parameters will accept the
  following arguments: \code{main} an overall title for the
  plot.  \code{xlab, ylab} a title for the x and yaxis.
  \code{xlim, ylim} numeric vectors of length 2, giving the
  x and y coordinates ranges.  \code{pch} either an integer
  specifying a symbol or a single character to be used as
  the default in plotting points.}
}
\description{
  \code{plot.iNEXT} the \code{\link{plot}} method for
  \code{\link{iNEXT}} Object
}
\examples{
data(spider)
x <- iNEXT.Ind(spider$Girdled, q=0)
plot(x)
data(ant)
y <- iNEXT.Sam(ant$h500m, q=1, t=round(seq(10, 500, length.out=20)), se=FALSE)
plot(y)
}
\seealso{
  \code{\link{lines.iNEXT}}
}

