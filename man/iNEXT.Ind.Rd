\name{iNEXT.Ind}
\alias{iNEXT.Ind}
\title{iNterpolation and EXTrapolation of abundance-based Hill number}
\usage{
  iNEXT.Ind(Spec, q = 0, m = NULL,
    endpoint = 2 * sum(Spec), Knots = 40, se = TRUE,
    nboot = 200)
}
\arguments{
  \item{Spec}{a vector of species abundances}

  \item{q}{a numeric value, the order of Hill number}

  \item{m}{a integer vector of rarefaction/extrapolation
  sample size, default is NULL. If m is not be specified,
  then the program will compute sample units due to
  endpoint and knots.}

  \item{endpoint}{a integer of sample size that is the
  endpoint for rarefaction/extrapolation. Default is double
  the original sample size.}

  \item{Knots}{a number of knots of computation, default is
  40}

  \item{se}{calculate bootstrap standard error and 95%
  confidence interval; default is TRUE}

  \item{nboot}{the number of bootstrap resampling times,
  default is 200}
}
\value{
  a list of interpolation and extrapolation Hill number
  with specific order q (qD) and sample coverage (SC)
}
\description{
  \code{iNEXT.Ind} Estimation of interpolation and
  extrapolation of abundance-based Hill number with order q
}
\examples{
data(spider)
# q = 0 with specific endpoint
iNEXT.Ind(spider$Girdled, q=0, endpoint=500)
# q = 1 with specific sample size m and don't calculate standard error
iNEXT.Ind(spider$Girdled, q=1, m=c(1, 10, 20, 50, 100, 200, 400, 600), se=FALSE)
}
\seealso{
  \code{\link{iNEXT.Sam}}
}

