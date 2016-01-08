\name{plant}
\alias{plant}
\docType{data}
\title{
Three coastal dune plant communities
}
\description{
The datset contains 272 random vegetation plots sampled in the
period 2002–2009 in the three most typical foredune habitats:
(i) upper beach and embryo dunes (\code{Embryo}), (ii) mobile dunes and (\code{Mobile}) (iii) transition dunes (\code{Transition}). For further details see Carboni, Santoro & Acosta (2010, 2011).
}
\usage{data("plant")}
\format{
  A list 3: \cr
  $Embryo is a data frame with 46 species (rows) on the following 70 sites (columns). \cr
  $Mobile is a data frame with 46 species (rows) on the following 131 sites (columns). \cr
  $Transition is a data frame with 46 species (rows) on the following 71 sites (columns). 
}
\details{
Random sampling was performed using georeferenced points extending over 250 km along the Tyrrhenian coast (Lazio Region; from 42º23'N, 11º39'E to 41º11'N, 13º20'E) and distributed in several study sites including most of the best conserved remnant dune systems of the region (for further details see Carboni, Santoro & Acosta 2010, 2011). At each random point, all vascular plants sampled in a square 2 x 2 mplot were listed together with a visual estimate of each species using a 10\%-interval rank scale. For each community type, all plant species present in at least 20\% of plots (for a total of 47 species) were described by a set of 16 functional traits related to different plant functions.
}
\source{
}
\references{
Carboni, M., Santoro, R. & Acosta, A.T.R. (2010) Are some communities of the coastal dune zonation more susceptible to alien plant invasion? Journal of Plant Ecology, 3, 139–147. \cr

Carboni, M., Santoro, R. & Acosta, A.T. (2011) Dealing with scarse data to understand how environmental gradients and propagule pressure shape finescale distribution patterns on coastal dunes. Journal of Vegetation Science, 22, 751–765. \cr

Ricotta, C., Pavoine, S., Bacaro, G. & Acosta, A.T.R. (2012) Functional rarefaction for species abundance data. Methods in Ecology and Evolution, 3, 519-525.
}
\examples{
data(plant)
out <- iNEXT(plant, datatype="abundance", rowsum=FALSE)
ggiNEXT(out)
}
\keyword{datasets}
