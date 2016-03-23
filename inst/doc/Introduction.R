## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "", 
                      fig.retina=2,
                      fig.align='center',
                      fig.width = 8, fig.height = 6,
                      out.width = "100%")
library(iNEXT)
library(ggplot2)
data("spider")
data("ant")
data("plant")

## ----eval=FALSE----------------------------------------------------------
#  ## install iNEXT package from CRAN
#  install.packages("iNEXT")
#  
#  ## install the latest version from github
#  install.packages('devtools')
#  library(devtools)
#  install_github('JohnsonHsieh/iNEXT')
#  
#  ## import packages
#  library(iNEXT)
#  library(ggplot2)

## ----eval=FALSE----------------------------------------------------------
#  iNEXT(x, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, nboot=50)

## ----eval=FALSE----------------------------------------------------------
#  data(spider)
#  str(spider)
#  iNEXT(spider, q=0, datatype="abundance")

## ----eval=FALSE----------------------------------------------------------
#  # set a series of sample sizes (m) for R/E computation
#  m <- c(1, 5, 20, 50, 100, 200, 400)
#  iNEXT(spider, q=0, datatype="abundance", size=m)

## ----eval=FALSE----------------------------------------------------------
#  out <- iNEXT(spider, q=c(0,1,2), datatype="abundance", size=m)

## ----eval=FALSE----------------------------------------------------------
#  ggiNEXT(x, type=1, se=TRUE, facet.var="none", color.var="site", grey=FALSE)

## ----eval=FALSE----------------------------------------------------------
#  out <- iNEXT(spider, q=c(0, 1, 2), datatype="abundance", endpoint=500)
#  # Sample‐size‐based R/E curves, separating by "site""
#  ggiNEXT(out, type=1, facet.var="site")
#  ## Not run:
#  # Sample‐size‐based R/E curves, separating by "order"
#  ggiNEXT(out, type=1, facet.var="order")
#  # display black‐white theme
#  ggiNEXT(out, type=1, facet.var="order", grey=TRUE)
#  ## End(Not run)

## ----eval=FALSE----------------------------------------------------------
#  # Sample‐size‐based R/E curves, separating by "site""
#  ggiNEXT(out, type=1, facet.var="site")

## ----echo=FALSE----------------------------------------------------------
out <- iNEXT(spider, q=c(0, 1, 2), datatype="abundance", endpoint=500)
ggiNEXT(out, type=1, facet.var="site")

## ------------------------------------------------------------------------
ggiNEXT(out, type=1, facet.var="order", color.var="site")

## ------------------------------------------------------------------------
ggiNEXT(out, type=2, facet.var="none", color.var="site")

## ------------------------------------------------------------------------
ggiNEXT(out, type=3, facet.var="site")

## ------------------------------------------------------------------------
ggiNEXT(out, type=3, facet.var="order", color.var="site")

## ------------------------------------------------------------------------
data(ant)
str(ant)

## ------------------------------------------------------------------------
t <- seq(1, 700, by=10)
out.inc <- iNEXT(ant, q=0, datatype="incidence_freq", size=t)

# Sample‐size‐based R/E curves
ggiNEXT(out.inc, type=1, color.var="site") + 
  theme_bw(base_size = 18) + 
  theme(legend.position="none")

## ------------------------------------------------------------------------
# Sample completeness curves
ggiNEXT(out.inc, type=2, color.var="site") +
  ylim(c(0.9,1)) +
  theme_bw(base_size = 18) + 
  theme(legend.position="none")

## ------------------------------------------------------------------------
# Coverage‐based R/E curves
ggiNEXT(out.inc, type=3, color.var ="site") + 
  xlim(c(0.9,1)) +
  theme_bw(base_size = 18) +
  theme(legend.position="bottom",
        legend.title=element_blank())

## ----eval=FALSE----------------------------------------------------------
#  estimateD(x, datatype="abundance", base="size", level=NULL)

## ------------------------------------------------------------------------
estimateD(ant, datatype="incidence_freq", 
          base="coverage", level=0.985)

## ------------------------------------------------------------------------
data(plant)
str(plant)

## ------------------------------------------------------------------------
out.raw <- iNEXT(plant, datatype="incidence_raw", endpoint=150)
ggiNEXT(out.raw)

## ------------------------------------------------------------------------
ggiNEXT(out, type=3, facet.var="site") + 
  theme(legend.position="none")

## ------------------------------------------------------------------------
ggiNEXT(out, type=1, facet.var="site") + 
  theme_bw(base_size = 18) +
  theme(legend.position="right")

## ------------------------------------------------------------------------
ggiNEXT(out, type=1, facet.var="order", grey=TRUE)

## ------------------------------------------------------------------------
ggiNEXT(out, type=1, facet.var="order") + 
  facet_wrap(~order, scales="free")

## ------------------------------------------------------------------------
ggiNEXT(out, type=1, facet.var="site") +
  scale_shape_manual(values=c(19,19,19))

