\name{Chat.Sam}
\alias{Chat.Sam}
\title{Incidence-based sample coverage}
\usage{
  Chat.Sam(x, t)
}
\arguments{
  \item{x}{a vector of species incidence-based frequency,
  the first entry is the total number of sampling units,
  followed by the speceis incidences abundances.}

  \item{t}{a integer vector of rarefaction/extrapolation
  sample size}
}
\value{
  a vector of estimated sample coverage function
}
\description{
  \code{Chat.Sam} Estimation of incidence-based sample
  coverage function
}

