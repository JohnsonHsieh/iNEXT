#'---
#'title: "Analysis species abundance data by iNEXT"
#'author: "T.C. Hsieh"
#'date: "Wednesday, August 27, 2014"
#'output: pdf_document
#'---

#'This is an simple guide for iNEXT. First install iNEXT R package from github repository (iNEXT is prepareing to submit to CRAN) by the following commends:

# install.packages('devtools') # Tools to install R packages from github github repositories
devtools::install_github('iNEXT', 'JohnsonHsieh') # install iNEXT version 2.0 by devtools

#' Next, import iNEXT package and dataset

library(iNEXT)
spec <- cbind(High =c(7,0,9,19,84,159),
              Intermediate =c(0,0,45,41,146,63),
              Low =c(15,3,3,10,126,4))
spec

#'Third, calculus species accumulation for each sites
#' Compute species accumulations
out <- iNEXT(x=spec, datatype="abundance", endpoint=200, knots=100, nboot=100)
#' Data visualization
#+ fig.width=8, fig.height=7
ggiNEXT(out, type=1, color.var = "site")
ggiNEXT(out, type=2, color.var = "site")
ggiNEXT(out, type=3, color.var = "site")
ggiNEXT(out, type=1, facet.var = "site", color.var = "site")

#' Further, you can compute species diversity accumulations (e.g. Shannon entropy and Simpson index), see more introduction in [Chao et al. (2014)](http://www.esajournals.org/doi/abs/10.1890/13-0133.1)
#' Compute species accumulations
out <- iNEXT(x=spec, q=c(0,1,2), datatype="abundance", endpoint=200, knots=100, nboot=100)
#' Data visualization
#+ fig.width=8, fig.height=7
ggiNEXT(out, type=1, facet.var = "site", color.var = "order")
ggiNEXT(out, type=2, color.var = "site")
ggiNEXT(out, type=3, facet.var = "order", color.var = "both")


#' ## How to cite this package

#' Use the following command
#+ comment=""
citation("iNEXT")
