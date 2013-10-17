\name{EstDis}
\alias{EstDis}
\title{Estimation of the rank of species relative abundance distribution or detection probability}
\usage{
  EstDis(x, datatype = c("abundance", "incidence"))
}
\arguments{
  \item{x}{a vector of species abundance or incidence
  frequency. If \code{datatype = "incidence"}, then the
  input format of first entry should be total number of
  sampling units, and followed by species incidence
  frequency.}

  \item{datatype}{the data type of input data. That is
  individual-based abundance data (\code{datatype =
  "abundance"}) or sample-based incidence data
  (\code{datatype = "incidence"}).}
}
\value{
  a vector of the rank of estimated relative abundance
  distribution or detection probability
}
\description{
  \code{EstDis} Estimation of the rank of species relative
  abundance distribution or detection probability to obtain
  bootstrap s.e.
}
\examples{
data(spider)
EstDis(spider$Girdled, datatype="incidence")
data(ant)
EstDis(ant$h50m, datatype="incidence")
}

