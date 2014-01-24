\name{iNEXT}
\alias{iNEXT}
\title{iNterpolation and EXTrapolation of Hill number}
\usage{
  iNEXT(x, q = 0, datatype = c("abundance", "incidence"),
    size = NULL, endpoint = NULL, knots = 40, se = TRUE,
    nboot = 50)
}
\arguments{
  \item{x}{a vector of species abundance or incidence
  frequency. If \code{datatype = "incidence"}, then the
  input format of first entry should be total number of
  sampling units, and followed by species incidence
  frequency.}

  \item{q}{a numeric value, the order of Hill number .}

  \item{datatype}{the data type of input data. That is
  individual-based abundance data (\code{datatype =
  "abundance"}) or sample-based incidence data
  (\code{datatype = "incidence"}).}

  \item{size}{an integer vector of
  rarefaction/extrapolation sample size (number of
  individuals or sampling units), default is NULL. If
  \code{size} is not be specified, it would compute
  rarefaction/extrapolation by endpoint and knots.}

  \item{endpoint}{a integer of sample size that is the
  endpoint for rarefaction/extrapolation. Default is double
  reference sample size.}

  \item{knots}{a number of knots (say K, default is 40)
  specifying separate sample size between 1 and
  \code{endpoint}. If \code{endpoint} is smaller than
  reference sample size, then \code{iNEXT()} compute
  rarefaction part only and divided by approximately equall
  spaced K. If \code{endpoint} is larger than reference
  sample size, then \code{iNEXT()} will compute
  approximately K/2 equally spaced for rarefaction part and
  divided extrapolation part as approximately K/2 equally
  spaced between reference sample size and
  \code{endpoint}.}

  \item{se}{calculate bootstrap standard error and show 95%
  confidence interval; default is TRUE.}

  \item{nboot}{the number of bootstrap resampling times,
  default is 50.}
}
\value{
  a list of interpolation and extrapolation Hill number
  with specific order q (qD) and sample coverage (SC)
}
\description{
  \code{iNEXT} Estimation of interpolation and
  extrapolation of Hill number with order q
}
\examples{
data(spider)
iNEXT(spider$Girdled, q=0, datatype="abundance")
data(ant)
iNEXT(ant$h500m, q=1, datatype="incidence", size=round(seq(10, 500, length.out=20)), se=FALSE)
}

