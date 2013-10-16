\name{lines.iNEXT}
\alias{lines.iNEXT}
\title{Lines Method for an iNEXT Object}
\usage{
  lines.iNEXT(x, style = "N2D", col = 1, ...)
}
\arguments{
  \item{x}{a \code{\link{iNEXT}} objext computed by
  \code{\link{iNEXT.Ind}} or \code{\link{iNEXT.Sam}}}

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
data(ant)
x50 <- iNEXT.Sam(ant$h50m, q=0)
x500 <- iNEXT.Sam(ant$h500m, q=0)
plot(x50, ylim=c(1, 300))
lines(x500, col=2)
}
\seealso{
  \code{\link{lines.iNEXT}}
}

