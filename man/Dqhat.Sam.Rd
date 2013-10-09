\name{Dqhat.Sam}
\alias{Dqhat.Sam}
\title{iNterpolation and EXTrapolation of incidence-based Hill number}
\usage{
  Dqhat.Sam(y, q, t)
}
\arguments{
  \item{y}{a vector of species incidence-based frequency,
  the first entry is the total number of sampling units,
  followed by the speceis incidences abundances.}

  \item{q}{a numerical value of the order of Hill number}

  \item{t}{a integer vector of rarefaction/extrapolation
  sample size}
}
\value{
  a vector of estimated interpolation and extrapolation
  function of Hill number with order q
}
\description{
  \code{Dqhat.Sam} Estimation of interpolation and
  extrapolation of incidence-based Hill number
}

