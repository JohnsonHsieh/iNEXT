## ----include=FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "", 
                      fig.retina=2,
                      fig.align='center',
                      fig.width = 8, fig.height = 6,
                      out.width = "100%")

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

