\name{ant}
\alias{ant}
\docType{data}
\title{
Ant data (datatype = "incidence_freq")
}
\description{
Ant species incidence frequencies for samples from five elevations/assemblages in northeastern Costa Rica (Longino and Colwell 2011). The number of sampling units (1m x 1m forest floor plot) for the 5 assemblages are respectively 599, 230, 150, 200 and 200. The number of observed species for the 5 assemblages are respectively 227, 241, 122, 56 and 14.
}
\usage{data(ant)}
\format{
  The input format for each site is a list of incidence frequencies. For incidence data, the first entry must be the total number of sampling units, followed by the species incidence frequencies as shown below: \cr
A list of 5 vectors \cr
 $ h50m  : num [1:228] 599 1 1 1 1 1 1 1 1 1 ... \cr
 $ h500m : num [1:242] 230 1 1 1 1 1 1 1 1 1 ... \cr
 $ h1070m: num [1:123] 150 1 1 1 1 1 1 1 1 1 ... \cr
 $ h1500m: num [1:57] 200 1 1 1 1 1 1 1 1 1 ... \cr
 $ h2000m: num [1:15] 200 1 2 2 3 4 8 8 13 15 ... \cr
 }
\references{
Longino, J.T. & Colwell, R.K. (2011) Density compensation, species composition, and richness of ants on a neotropical elevational gradient. Ecosphere, 2, art29.
}

\keyword{datasets}
