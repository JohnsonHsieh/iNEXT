\name{lines.iNEXT}
\alias{lines.iNEXT}
\title{Lines Method for an iNEXT Object}
\usage{
  lines.iNEXT(x, style = "N2D", col = 1, ...)
}
\arguments{
  \item{x}{a \code{\link{iNEXT}} objext computed by
  \code{\link{iNEXT}}}

  \item{style}{three different plotting style = c("N2D",
  "N2SC", "SC2D"),}

  \item{col}{a specification for the default plotting
  color, see \code{\link{par}} for detail.}

  \item{...}{further plotting parameters (see
  \code{\link{par}}) may also be accepted.}
}
\description{
  \code{lines.iNEXT} the \code{\link{lines}} method for
  \code{\link{iNEXT}} Object
}
\examples{
data(spider)
x1 <- iNEXT(spider$Girdled, q=0, datatype="abundance")
x2 <- iNEXT(spider$Logged, q=0, datatype="abundance")
plot(x1, ylim=c(1, 300))
lines(x2, col=2)
}
\seealso{
  \code{\link{lines.iNEXT}}
}

