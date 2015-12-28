\name{BCI.inc}
\alias{BCI.inc}
\docType{data}
\title{
Barro Colorado Island Tree Counts
}
\description{
Tree counts in 1-hectare plots in the Barro Colorado Island and associated site information.
}
\usage{data("BCI.inc")}
\format{
  A data frame with 299 species (rows) on the following 200 sites (columns).
}
\details{
This data frame is a species by sites presence-absence matrix from the 50-ha Barro Colorado Island (BCI), Panama, 1985 tree census. Data give the number of trees and shrubs at least 1 cm in diameter at breast height.
}
\source{
http://ctfs.si.edu/
}
\references{
Hubbell, S.P., Condit, R., and Foster, R.B. 2005. Barro Colorado Forest Census Plot Data. URL http://ctfs.si.edu/webatlas/datasets/bci. 
}
\examples{
data(BCI.inc)
dat <- as.incfreq(BCI.inc)
DataInfo(dat, datatype = "incidence")
}
\keyword{datasets}
