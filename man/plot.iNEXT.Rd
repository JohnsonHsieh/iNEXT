\name{plot.iNEXT}
\alias{plot.iNEXT}
\title{Plot Method for an iNEXT Object}
\usage{
  plot.iNEXT(x, style = "N2D", col = 1, ...)
}
\arguments{
  \item{x}{a \code{\link{iNEXT}} objext computed by
  \code{\link{iNEXT}}}

  \item{style}{three different plotting style = c("N2D",
  "N2SC", "SC2D"); N2D means to plot number of individuals
  or number of samples to diversity; N2SC means to plot
  number of individuals or number of samples to sample
  coverage; SC2D means to plot sample coverage to
  diversity.}

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
x <- iNEXT(spider$Girdled, q=0, datatype="abundance")
plot(x)
data(ant)
y <- iNEXT(ant$h500m, q=1, datatype="incidence", size=round(seq(10, 500, length.out=20)), se=FALSE)
plot(y)
}
\seealso{
  \code{\link{lines.iNEXT}}
}

