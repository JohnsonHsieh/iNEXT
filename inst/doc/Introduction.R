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
#  iNEXT(x, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)

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
#  data(bird)
#  str(bird) # 41 obs. of 2 variables
#  iNEXT(bird, q=0, datatype="abundance")

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

## ------------------------------------------------------------------------
library(iNEXT)
library(ggplot2)
library(gridExtra)
library(grid)
data("spider")
out <- iNEXT(spider, q=0, datatype="abundance")
g <- ggiNEXT(out, type=1, color.var = "site")
g

## ------------------------------------------------------------------------
g1 <- g + scale_shape_manual(values=c(11, 12)) + 
          scale_linetype_manual(values=c(1,2))
g2 <- g + scale_colour_manual(values=c("red", "blue")) +
          scale_fill_manual(values=c("red", "blue"))

# Draw multiple graphical objec on a page
# library(gridExtra)
grid.arrange(g1, g2, ncol=2)

## ------------------------------------------------------------------------
# point is drawn on the 1st layer, default size is 5
gb3 <- ggplot_build(g)
gb3$data[[1]]$size <- 10
gt3 <- ggplot_gtable(gb3)

# use grid.draw to draw the graphical object
# library(grid)
# grid.draw(gt3)

## ------------------------------------------------------------------------
# line is drawn on the 2nd layer, default size is 1.5
gb4 <- ggplot_build(g)
gb4$data[[2]]$size <- 3
gt4 <- ggplot_gtable(gb4)
# grid.draw(gt4)

## ------------------------------------------------------------------------
grid.arrange(gt3, gt4, ncol=2)

## ------------------------------------------------------------------------
g5 <- g + theme_bw() + theme(legend.position = "bottom")
g6 <- g + theme_classic() + theme(legend.position = "bottom")
grid.arrange(g5, g6, ncol=2)

## ------------------------------------------------------------------------
library(ggthemes)
g7 <- g + theme_hc(bgcolor = "darkunica") +
          scale_colour_hc("darkunica")

g8 <- g + theme_economist() + scale_colour_economist()
grid.arrange(g7, g8, ncol=2)

## ------------------------------------------------------------------------
g9 <- g + theme_bw(base_size = 18) +
      scale_fill_grey(start = 0, end = .4) +
      scale_colour_grey(start = .2, end = .2) +
      theme(legend.position="bottom",
            legend.title=element_blank())

g10 <- g + theme_tufte(base_size = 12) +       
    scale_fill_grey(start = 0, end = .4) +
    scale_colour_grey(start = .2, end = .2) +
    theme(legend.position="bottom",
          legend.title=element_blank())
grid.arrange(g9, g10, ncol=2)

## ------------------------------------------------------------------------
df <- fortify(out, type=1)
head(df)

df.point <- df[which(df$method=="observed"),]
df.line <- df[which(df$method!="observed"),]
df.line$method <- factor(df.line$method, 
                         c("interpolated", "extrapolated"),
                         c("interpolation", "extrapolation"))
 
ggplot(df, aes(x=x, y=y, colour=site)) + 
  geom_point(aes(shape=site), size=5, data=df.point) +
  geom_line(aes(linetype=method), lwd=1.5, data=df.line) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=site, colour=NULL), alpha=0.2) +
  labs(x="Number of individuals", y="Species diversity") +
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=18)) 

