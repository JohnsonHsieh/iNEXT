\name{bird}
\alias{bird}
\docType{data}
\title{Bird data (datatype = "abundance")}
\description{
This data set includes the abundances of 41 bird species collected in two sites (the North and South sites) at the Barrington Tops National Park, Australia (Chao et al. 2015)
}
\usage{data(bird)}
\format{
  a data.frame with 41 species (rows) and two sites (columns).
}
\source{
Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015) Rarefaction and extrapolation of phylogenetic diversity. Methods in Ecology and Evolution, 6, 380-388. 
}
\examples{
data(bird)
\dontrun{
out <- iNEXT(bird, datatype="abundance")
ggiNEXT(out)
}
}
\keyword{datasets}
