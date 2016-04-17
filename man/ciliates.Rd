\name{ciliates}
\alias{ciliates}
\docType{data}
\title{
Three soil samples
}
\description{
A total of 51 soil samples were taken from three areas ( \code{EtoshaPan}, \code{CentralNamibDesert}, \code{SouthernNamibDesert}) of Namibia. In area \code{EtoshaPan}, there were 19 soil samples and the number of observed species was 216. In area \code{CentralNamibDesert}, there were 17 soil samples and the number of observed species was 130. In area \code{SouthernNamibDesert}, there were 15 soil samples and the number of observed species was 150.
}
\usage{data("ciliates")}
\format{
  A list 3: \cr
  $EtoshaPan is a matrix with 365 species (rows) on the following 19 sites (columns). \cr
  $CentralNamibDesert is a matrix with 365 species (rows) on the following 17 sites (columns). \cr
  $SouthernNamibDesert is a matrix with 365 species (rows) on the following 15 sites (columns). 
}

\references{
Foissner, W., Agatha, S., & Berger, H. (2002) Soil Ciliates (Protozoa, Ciliophora) From Namibia (Southwest Africa), With Emphasis on Two Contrasting Environments, the Etosha Region and the Namib Desert. Denisia, 5, 1-1459. 

}
\examples{
data(ciliates)
out <- iNEXT(ciliates, datatype="incidence_raw")
ggiNEXT(out)
}
\keyword{datasets}
