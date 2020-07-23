
#
#
###############################################
#' ggplot2 extension for an iNEXT object
#' 
#' \code{ggiNEXT}: the \code{\link[ggplot2]{ggplot}} extension for \code{\link{iNEXT}} Object to plot sample-size- and coverage-based rarefaction/extrapolation curves along with a bridging sample completeness curve
#' @param x an \code{iNEXT} object computed by \code{\link{iNEXT}}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); unconditional coverage-based rarefaction/extrapolation curve (\code{type = 3}); conditional coverage-based rarefaction/extrapolation curve (\code{type = 4}).                
#' @param se a logical variable to display confidence interval around the estimated sampling curve.
#' @param facet.var create a separate plot for each value of a specified variable: 
#'  no separation \cr (\code{facet.var="None"}); 
#'  a separate plot for each diversity order (\code{facet.var="Order.q"}); 
#'  a separate plot for each assemblage (\code{facet.var="Assemblage"}); 
#'  a separate plot for each combination of order x assemblage (\code{facet.var="Both"}).              
#' @param color.var create curves in different colors for values of a specified variable:
#'  all curves are in the same color (\code{color.var="None"}); 
#'  use different colors for diversity orders (\code{color.var="Order.q"}); 
#'  use different colors for sites (\code{color.var="Assemblage"}); 
#'  use different colors for combinations of order x assemblage (\code{color.var="Both"}).  
#' @param grey a logical variable to display grey and white ggplot2 theme. 
#' @param ... other arguments passed on to methods. Not currently used.
#' @return a ggplot2 object
#' @examples
#' data(spider)
#' # single-assemblage abundance data
#' out1 <- iNEXT(spider$Girdled, q=0, datatype="abundance")
#' ggiNEXT(x=out1, type=1)
#' ggiNEXT(x=out1, type=2)
#' ggiNEXT(x=out1, type=3)
#' 
#'\dontrun{
#' # single-assemblage incidence data with three orders q
#' data(ant)
#' size <- round(seq(10, 500, length.out=20))
#' y <- iNEXT(ant$h500m, q=c(0,1,2), datatype="incidence_freq", size=size, se=FALSE)
#' ggiNEXT(y, se=FALSE, color.var="Order.q")
#' 
#' # multiple-assemblage abundance data with three orders q
#' z <- iNEXT(spider, q=c(0,1,2), datatype="abundance")
#' ggiNEXT(z, facet.var="Assemblage", color.var="Order.q")
#' ggiNEXT(z, facet.var="Both", color.var="Both")
#'}
#' @export
#' 
ggiNEXT <- function(x, type=1, se=TRUE, facet.var="None", color.var="Assemblage", grey=FALSE){  
  UseMethod("ggiNEXT", x)
}

#' @export
#' @rdname ggiNEXT
ggiNEXT.iNEXT <- function(x, type=1, se=TRUE, facet.var="None", color.var="Assemblage", grey=FALSE){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#330066", "#CC79A7",  "#0072B2", "#D55E00"))
  TYPE <-  c(1, 2, 3, 4)
  SPLIT <- c("None", "Order.q", "Assemblage", "Both")
  if(is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1)
    stop("invalid plot type")
  if(is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var, SPLIT) == -1)
    stop("invalid facet variable")
  if(is.na(pmatch(color.var, SPLIT)) | pmatch(color.var, SPLIT) == -1)
    stop("invalid color variable")
  
  type <- pmatch(type, 1:4)
  facet.var <- match.arg(facet.var, SPLIT)
  color.var <- match.arg(color.var, SPLIT)
  if(facet.var=="Order.q") color.var <- "Assemblage"
  if(facet.var=="Assemblage") color.var <- "Order.q"
  
  options(warn = -1)
  z <- fortify(x, type=type)
  options(warn = 0)
  if(!('y.lwr' %in% names(z))) { se <- FALSE }
  datatype <- unique(z$datatype)
  if(color.var=="None"){
    if(levels(factor(z$Order.q))>1 & length(unique(z$Assemblage))>1){
      warning("invalid color.var setting, the iNEXT object consists multiple assemblages and orders, change setting as Both")
      color.var <- "Both"
      z$col <- z$shape <- paste(z$Assemblage, z$Order.q, sep="-")
      
    }else if(length(unique(z$Assemblage))>1){
      warning("invalid color.var setting, the iNEXT object consists multiple assemblages, change setting as Assemblage")
      color.var <- "Assemblage"
      z$col <- z$shape <- z$Assemblage
    }else if(levels(factor(z$Order.q))>1){
      warning("invalid color.var setting, the iNEXT object consists multiple orders, change setting as Order.q")
      color.var <- "Order.q"
      z$col <- z$shape <- factor(z$Order.q)
    }else{
      z$col <- z$shape <- rep(1, nrow(z))
    }
  }else if(color.var=="Order.q"){     
    z$col <- z$shape <- factor(z$Order.q)
  }else if(color.var=="Assemblage"){
    if(length(unique(z$Assemblage))==1){
      warning("invalid color.var setting, the iNEXT object do not consist multiple assemblages, change setting as Order.q")
      z$col <- z$shape <- factor(z$Order.q)
    }
    z$col <- z$shape <- z$Assemblage
  }else if(color.var=="Both"){
    if(length(unique(z$Assemblage))==1){
      warning("invalid color.var setting, the iNEXT object do not consist multiple assemblages, change setting as Order.q")
      z$col <- z$shape <- factor(z$Order.q)
    }
    z$col <- z$shape <- paste(z$Assemblage, z$Order.q, sep="-")
  }
  zz=z
  z$Method[z$Method=="Observed"]="Rarefaction"
  z$lty <- factor(z$Method, levels = c("Rarefaction", "Extrapolation"))
  z$col <- factor(z$col)
  data.sub <- zz[which(zz$Method=="Observed"),]
  
  g <- ggplot(z, aes_string(x="x", y="y", colour="col")) + 
    geom_point(aes_string(shape="shape"), size=5, data=data.sub)+
    scale_colour_manual(values=cbPalette)
  
  
  g <- g + geom_line(aes_string(linetype="lty"), lwd=1.5) +
    guides(linetype=guide_legend(title="Method"),
           colour=guide_legend(title="Guides"), 
           fill=guide_legend(title="Guides"), 
           shape=guide_legend(title="Guides")) +
    theme(legend.position = "bottom", 
          legend.title=element_blank(),
          text=element_text(size=18),
          legend.key.width = unit(1.2,"cm")) 
  
  if(type==2L) {
    g <- g + labs(x="Number of sampling units", y="Sample coverage")
    if(datatype=="abundance") g <- g + labs(x="Number of individuals", y="Sample coverage")
  }else if(type==3L|type==4L) {
    g <- g + labs(x="Sample coverage", y="Species diversity")
  }else {
    g <- g + labs(x="Number of sampling units", y="Species diversity")
    if(datatype=="abundance") g <- g + labs(x="Number of individuals", y="Species diversity")
  }
  
  if(se)
    g <- g + geom_ribbon(aes_string(ymin="y.lwr", ymax="y.upr", fill="factor(col)", colour="NULL"), alpha=0.2)+
    scale_fill_manual(values=cbPalette)
  
  
  if(facet.var=="Order.q"){
    if(length(levels(factor(z$Order.q))) == 1 & type!=2){
      warning("invalid facet.var setting, the iNEXT object do not consist multiple orders.")      
    }else{
      odr_grp <- as_labeller(c(`0` = "q = 0", `1` = "q = 1",`2` = "q = 2")) 
      g <- g + facet_wrap(~Order.q, nrow=1, labeller = odr_grp)
      if(color.var=="Both"){
        g <- g + guides(colour=guide_legend(title="Guides", ncol=length(levels(factor(z$Order.q))), byrow=TRUE),
                        fill=guide_legend(title="Guides"))
      }
      if(type==2){
        g <- g + theme(strip.background = element_blank(),strip.text.x = element_blank())
          
      }
    }
  }
  
  if(facet.var=="Assemblage"){
    if(length(unique(z$Assemblage))==1) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple assemblages")
    }else{
      g <- g + facet_wrap(~Assemblage, nrow=1)
      if(color.var=="Both"){
        g <- g + guides(colour=guide_legend(title="Guides", nrow=length(levels(factor(z$Order.q)))),
                        fill=guide_legend(title="Guides"))
      }
    }
  }
  
  if(facet.var=="Both"){
    if(length(levels(factor(z$Order.q))) == 1 | length(unique(z$Assemblage))==1){
      warning("invalid facet.var setting, the iNEXT object do not consist multiple assemblages or orders.")
    }else{
      odr_grp <- as_labeller(c(`0` = "q = 0", `1` = "q = 1",`2` = "q = 2")) 
      g <- g + facet_wrap(Assemblage~Order.q,labeller = labeller(Order.q = odr_grp)) 
      if(color.var=="both"){
        g <- g +  guides(colour=guide_legend(title="Guides", nrow=length(levels(factor(z$Assemblage))), byrow=TRUE),
                         fill=guide_legend(title="Guides"))
      }
    }
  }
  
  if(grey){
    g <- g + theme_bw(base_size = 18) +
      scale_fill_grey(start = 0, end = .4) +
      scale_colour_grey(start = .2, end = .2) +
      guides(linetype=guide_legend(title="Method"), 
             colour=guide_legend(title="Guides"), 
             fill=guide_legend(title="Guides"), 
             shape=guide_legend(title="Guides")) +
      theme(legend.position="bottom",
            legend.title=element_blank())
  }
  g <- g + theme(legend.box = "vertical")
  return(g)
  
}

#' @export
#' @rdname ggiNEXT
ggiNEXT.default <- function(x, ...){
  stop(
    "iNEXT doesn't know how to deal with data of class ",
    paste(class(x), collapse = "/"),
    call. = FALSE
  )
}

#' Fortify method for classes from the iNEXT package.
#'
#' @name fortify.iNEXT
#' @param model \code{iNEXT} to convert into a dataframe.
#' @param data not used by this method
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); unconditional coverage-based rarefaction/extrapolation curve (\code{type = 3});
#' conditional coverage-based rarefaction/extrapolation curve (\code{type = 4}).                 
#' @param ... not used by this method
#' @export
#' @examples
#' data(spider)
#' # single-assemblage abundance data
#' out1 <- iNEXT(spider$Girdled, q=0, datatype="abundance")
#' ggplot2::fortify(out1, type=1)

fortify.iNEXT <- function(model, data = model$iNextEst, type = 1, ...) {
  datatype <- ifelse(names(model$DataInfo)[2]=="n","abundance","incidence")
  z <- data
  # if(class(z) == "list"){
  #   if(datatype=='abundance'){
  #     id_match <- match(c("Assemblage","m", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL", "SC"), names(z$coverage_based), nomatch = 0)  
  #   }else if (datatype=='incidence'){
  #     id_match <- match(c("Assemblage","t", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL", "SC"), names(z$coverage_based), nomatch = 0)  
  #   }
  #   z$coverage_based <- cbind(z$coverage_based[,id_match],SC.LCL=NA,SC.UCL=NA)
  #   z <- data.frame(do.call("rbind", z), base=rep(names(z), each = sapply(z, nrow)))
  #   rownames(z) <- NULL
  # }else{
  #   z$site <- ""
  # }
  
  # if(ncol(z)==6) {
  #   warning("invalid se setting, the iNEXT object do not consist confidence interval")
  #   se <- FALSE
  # }else if(ncol(z)>6) {
  #   se <- TRUE
  # }
  
  if(is.na(z$size_based$qD.LCL[1])) {
    warning("invalid se setting, the iNEXT object do not consist confidence interval")
    se <- FALSE
  }else{
    se <- TRUE
  }
  
  if(type==1L) {
    z <- z$size_based
    z$x <- z[,2]
    z$y <- z$qD
    z$datatype <- datatype
    z$plottype <- type
    if(se){
      z$y.lwr <- z$qD.LCL
      z$y.upr <- z$qD.UCL
      data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y","y.lwr","y.upr")]
    }else{
      data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y")]
    }
  }else if(type==2L){
    z <- z$size_based
    if(length(unique(z$Order.q))>1){
      z <- subset(z, Order.q==unique(z$Order.q)[1])
    }
    z$x <- z[,2]
    z$y <- z$SC
    z$datatype <- datatype
    z$plottype <- type
    if(se){
      z$y.lwr <- z$SC.LCL
      z$y.upr <- z$SC.UCL
      data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y","y.lwr","y.upr")]
    }else{
      data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y")]
    }
  }else if(type==3L){
    z <- z$coverage_based
    z$x <- z$SC
    z$y <- z$qD
    z$datatype <- datatype
    z$plottype <- type
    if(se){
      z$y.lwr <- z$qD.LCL
      z$y.upr <- z$qD.UCL
      data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y","y.lwr","y.upr")]
    }else{
      data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y")]
    }
  }else if(type==4L){
    z <- z$size_based
    z$x <- z$SC
    z$y <- z$qD
    z$datatype <- datatype
    z$plottype <- type
    if(se){
      z$y.lwr <- z$qD.LCL
      z$y.upr <- z$qD.UCL
      data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y","y.lwr","y.upr")]
    }else{
      data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y")]
    }
  }
  data
}

