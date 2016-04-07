###############################################
#' Plotting iNEXT object
#' 
#' \code{plot.iNEXT}: Plotting method for objects inheriting from class "iNEXT"
#' @param x an \code{iNEXT} object computed by \code{\link{iNEXT}}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).                 
#' @param se a logical variable to display confidence interval around the estimated sampling curve.
#' @param show.legend a logical variable to display legend.
#' @param show.main a logical variable to display title.
#' @param col a vector for plotting color
#' @param ... arguments to be passed to methods, such as graphical parameters (\code{\link{par}}).
#' @examples
#' data(spider)
#' # single-assemblage abundance data
#' out1 <- iNEXT(spider$Girdled, q=0, datatype="abundance")
#' plot(x=out1, type=1)
#' plot(x=out1, type=2)
#' plot(x=out1, type=3)
#' 

#' @export
plot.iNEXT <- function(x, type=1, se=TRUE, show.legend=TRUE, show.main=TRUE, col=NULL,...){
  
  if(class(x) != "iNEXT")
    stop("invalid object class")
  TYPE <-  c(1, 2, 3)
  SPLIT <- c("none", "order", "site", "both")
  if(is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1)
    stop("invalid plot type")
  
  type <- pmatch(type, 1:3)
  
  y <- method <- site <- shape <- y.lwr <- y.upr <- NULL
  site <<- NULL
  
  z <- x$iNextEst
  if(class(z) == "list"){
    z <- data.frame(do.call("rbind", z), site=rep(names(z), sapply(z, nrow)))
    rownames(z) <- NULL
  }else{
    z$site <- ""
    z$site <- factor(z$site)
  }
  
  if("qD.LCL" %in% names(z) == FALSE & se) {
    warning("invalid se setting, the iNEXT object do not consist confidence interval")
    se <- FALSE
  }else if("qD.LCL" %in% names(z) & se) {
    se <- TRUE
  }else{
    se <- FALSE
  }
  
  if(type==1L) {
    z$x <- z[,1]
    z$y <- z$qD
    if(!is.null(xlab)) xlab <- ifelse(names(x$DataInfo)[2]=="n", "Number of individuals", "Number of sampling units")
    if(!is.null(ylab)) ylab <- "Species diversity"
    if(se){
      z$y.lwr <- z$qD.LCL
      z$y.upr <- z$qD.UCL
    }
  }else if(type==2L){
    if(length(unique(z$order))>1){
      z <- subset(z, order==unique(z$order)[1])
    }
    z$x <- z[,1]
    z$y <- z$SC
    if(!is.null(xlab)) xlab <- ifelse(names(x$DataInfo)[2]=="n", "Number of individuals", "Number of sampling units")
    if(!is.null(ylab)) ylab <- "Sample coverage"
    if(se){
      z$y.lwr <- z$SC.LCL
      z$y.upr <- z$SC.UCL
    }
  }else if(type==3L){
    z$x <- z$SC
    z$y <- z$qD
    if(!is.null(xlab)) xlab <- "Sample coverage"
    if(!is.null(ylab)) ylab <- "Species diversity"
    if(se){
      z$y.lwr <- z$qD.LCL
      z$y.upr <- z$qD.UCL
    }
  }
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  
  conf.reg=function(x,LCL,UCL,...) {
    x.sort <- order(x)
    x <- x[x.sort]
    LCL <- LCL[x.sort]
    UCL <- UCL[x.sort]
    polygon(c(x,rev(x)),c(LCL,rev(UCL)), ...)
  }
  
  SITE <- levels(z$site)
  ORDER <- unique(z$order)
  
  if(is.null(col)){
    col <- gg_color_hue(length(SITE))
  }else{
    col <- rep(col,length(SITE))[1:length(SITE)]
  }
  pch <- (16+1:length(SITE))%%25
  
  for(j in 1:length(ORDER)){
    if(se==TRUE){
      tmp.sub <- subset(z, order==ORDER[j])
      tmp.j <- data.frame(site=tmp.sub$site, order=tmp.sub$order,
                          method=tmp.sub$method, 
                          x=tmp.sub$x, y=tmp.sub$y,
                          y.lwr=tmp.sub$y.lwr, y.upr=tmp.sub$y.upr)
      
      plot(y.upr~x, data=tmp.j, type="n", xlab="", ylab="", ...)
    }else{
      tmp.sub <- subset(z, order==ORDER[j])
      
      tmp.j <- data.frame(site=tmp.sub$site, order=tmp.sub$order,
                          method=tmp.sub$method, 
                          x=tmp.sub$x, y=tmp.sub$y)
      
      plot(y~x, data=tmp.j, type="n", xlab="", ylab="", ...)
    }
    
    for(i in 1:length(SITE)){
      tmp <- subset(tmp.j, site==SITE[i])
      if(se==TRUE){
        conf.reg(x=tmp$x, LCL=tmp$y.lwr, UCL=tmp$y.upr, border=NA, col=adjustcolor(col[i], 0.25))
      }
      lines(y~x, data=subset(tmp, method=="interpolated"), lty=1, lwd=2, col=col[i])
      lines(y~x, data=subset(tmp, method=="extrapolated"), lty=2, lwd=2, col=col[i])
      points(y~x, data=subset(tmp, method=="observed"), pch=pch[i], cex=2, col=col[i])
      
    }
    if(show.legend==TRUE){
      if(type==3L){
        legend("topleft", legend=paste(SITE), col=col, lty=1, lwd=2, pch=pch, cex=1, bty="n")
      }else{
        legend("bottomright", legend=paste(SITE), col=col, lty=1, lwd=2, pch=pch, cex=1, bty="n")
      }
    }
    title(xlab=xlab, ylab=ylab)
    if(show.main==TRUE) title(main=paste("Order q =", ORDER[j]))
    par(ask=TRUE)
  }
  par(ask=FALSE)
}


#' Printing iNEXT object
#' 
#' \code{print.iNEXT}: Print method for objects inheriting from class "iNEXT"
#' @param x an \code{iNEXT} object computed by \code{\link{iNEXT}}.
#' @param ... additional arguments.
#' @export
print.iNEXT <- function(x, ...){
  site.n <- nrow(x$DataInfo)
  order.n <- ifelse(site.n > 1, 
                    paste(unique(x$iNextEst[[1]]$order), collapse = ", "),
                    paste(unique(x$iNextEst$order), collapse = ", "))
  cat("Compare ", site.n, " assemblages with Hill number order q = ", order.n,".\n", sep="")
  cat("$class: iNEXT\n\n")
  cat("$DataInfo: basic data information\n")
  print(x$DataInfo)
  cat("\n")
  cat("$iNextEst: diversity estimates with rarefied and extrapolated samples.\n")
  if(class(x$iNextEst)=="data.frame"){
    y <- x$iNextEst
    m <- quantile(y[,1], type = 1)
    res <- y[y[,1]%in%m,]
  }else{
    res <- lapply((x$iNextEst), function(y){
      m <- quantile(y[,1], type = 1)
      y[y[,1]%in%m,]
    })
  }
  print(res)
  cat("\n")
  cat("$AsyEst: asymptotic diversity estimates along with related statistics.\n")
  print(x$AsyEst)
  cat("\n")
  cat("NOTE: Only show five estimates, call iNEXT.objext$iNextEst. to show complete output.\n")
  return(invisible())
}
