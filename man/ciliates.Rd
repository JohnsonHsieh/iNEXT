\name{ciliates}
\alias{ciliates}
\docType{data}
\title{
Ciliates data (datatype = "incidence_raw") 
}
\description{
A total of 51 soil samples were taken from three areas ( \code{EtoshaPan}, \code{CentralNamibDesert}, \code{SouthernNamibDesert}) of Namibia. In area \code{EtoshaPan}, there were 19 soil samples and the number of observed species was 216. In area \code{CentralNamibDesert}, there were 17 soil samples and the number of observed species was 130. In area \code{SouthernNamibDesert}, there were 15 soil samples and the number of observed species was 150. The total number of species in the three areas was 365. The data include three lists of matrices; each matrix is a species by soil-sample matrix ("1" for a detection, and "0" otherwise).  
}
\usage{data("ciliates")}
\format{
  A list 3 matrices: \cr
  $EtoshaPan is a matrix with 365 species (rows) and 19 soil samples (columns). \cr
  $CentralNamibDesert is a matrix with 365 species (rows) and 17 soil samples (columns). \cr
  $SouthernNamibDesert is a matrix with 365 species (rows) and 15 soil samples (columns). 
}

\references{
Foissner, W., Agatha, S., & Berger, H. (2002) Soil Ciliates (Protozoa, Ciliophora) from Namibia (Southwest Africa), With Emphasis on Two Contrasting Environments, the Etosha Region and the Namib Desert. Denisia, 5, 1-1459. 

}
\examples{
data(ciliates)
\dontrun{
out <- iNEXT(ciliates, datatype = "incidence_raw")
ggiNEXT(out)
}
}
\keyword{datasets}
