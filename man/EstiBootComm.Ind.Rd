\name{EstiBootComm.Ind}
\alias{EstiBootComm.Ind}
\title{Estimation of species reletive abundance distribution}
\usage{
  EstiBootComm.Ind(Spec)
}
\arguments{
  \item{Spec}{a vector of species abundances}
}
\value{
  a vector of reltavie abundance
}
\description{
  \code{EstiBootComm.Ind} Estimation of species reletive
  abundance distribution to obtain bootstrap s.e.
}
\examples{
data(spider)
EstiBootComm.Ind(spider$Girdled)
}
\seealso{
  \code{\link{EstiBootComm.Sam}}
}

