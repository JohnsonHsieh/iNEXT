\name{bird}
\alias{bird}
\docType{data}
\title{Bird abundance-based data}
\description{
This data set describes 41 birds as reported by Jetz et al. (2012). It gives the two sites of species abundance to these 41 species in November 2012 at Barrington Tops National Park, Australia.
}
\usage{data(bird)}
\format{
  a data frame with 41 species and two sites: North and South sites.
}
\source{
Jetz, W., Thomas, G.H., Joy, J.B., Hartmann, K. & Mooers A.O. (2012). The global diversity of birds in space and time. Nature, 491, 444–448.
}
\examples{
library(ggplot2)
data(bird)
out <- iNEXT(bird, datatype="abundance")
ggiNEXT(out)
}
\keyword{datasets}
