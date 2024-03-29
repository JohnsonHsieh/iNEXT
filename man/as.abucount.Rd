% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/invChat.R
\name{as.abucount}
\alias{as.abucount}
\title{Transform abundance raw data to abundance row-sum counts (iNEXT input format)}
\usage{
as.abucount(x)
}
\arguments{
\item{x}{a \code{data.frame} or \code{matirx} (species by sites matrix).}
}
\value{
a \code{vector} of species abundance row-sum counts.
}
\description{
\code{as.abucount}: transform species abundance raw data (a species by sites matrix) to row-sum counts (iNEXT input format) as species abundances.
}
\examples{
data(ciliates)
lapply(ciliates, as.abucount)

}
