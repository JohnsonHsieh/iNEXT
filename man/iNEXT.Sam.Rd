\name{iNEXT.Sam}
\alias{iNEXT.Sam}
\title{iNterpolation and EXTrapolation of incidence-based Hill number}
\usage{
  iNEXT.Sam(Spec, t = NULL, q = 0,
    endpoint = 2 * max(Spec), Knots = 40, se = TRUE,
    nboot = 200)
}
\arguments{
  \item{Spec}{a vector of species incidence-based
  frequency, the first entry is the total number of
  sampling units, followed by the speceis incidences
  abundances.}

  \item{q}{a numeric value, the order of Hill number}

  \item{t}{a integer vector of rarefaction/extrapolation
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
  \code{iNEXT.Sam} Estimation of interpolation and
  extrapolation of incidence-based Hill number with order q
}
\examples{
data(ant)
# q = 0 with specific endpoint
iNEXT.Sam(ant$h50m, q=0, endpoint=100)
# q = 1 with specific sample size m and don't calculate standard error
iNEXT.Sam(ant$h500m, q=1, t=round(seq(10, 500, length.out=20)), se=FALSE)
}
\seealso{
  \code{\link{iNEXT.Ind}}
}

