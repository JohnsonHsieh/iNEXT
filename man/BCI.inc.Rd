\name{BCI.inc}
\alias{BCI.inc}
\docType{data}
\title{
Barro colorado island tree incidence data
}
\description{
The 50-hectare permanent tree plot was established in 1985 in the tropical moist forest of Barro Colorado Island (BCI) in Gatun Lake in central Panama. 
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
