\name{EstiBootComm.Sam}
\alias{EstiBootComm.Sam}
\title{Estimation of species detection distribution}
\usage{
  EstiBootComm.Sam(Spec)
}
\arguments{
  \item{Spec}{a vector of species incidence, the first
  entry is the total number of sampling units, followed by
  the speceis incidences abundances.}
}
\value{
  a vector of estimated detection probability
}
\description{
  \code{EstiBootComm.Sam} Estimation of species detection
  distribution to obtain bootstrap s.e.
}
\examples{
data(ant)
EstiBootComm.Sam(ant$h50m)
}
\seealso{
  \code{\link{EstiBootComm.Sam}}
}

