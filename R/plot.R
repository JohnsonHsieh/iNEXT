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
#' @importFrom grDevices adjustcolor
#' @importFrom grDevices hcl
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics points
#' @importFrom graphics polygon
#' @importFrom graphics title
#' @importFrom graphics plot
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
  
  if(!inherits(x, "iNEXT"))
    stop("invalid object class")
  TYPE <-  c(1, 2, 3)
  # SPLIT <- c("none", "order", "site", "both")
  if(is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1)
    stop("invalid plot type")
  
  type <- pmatch(type, 1:3)
  
  y <- method <- site <- shape <- y.lwr <- y.upr <- NULL
  site <<- NULL
  
  # z <- x$iNextEst
  # if(inherits(z, "list")){
  #   z <- data.frame(do.call("rbind", z), site=rep(names(z), sapply(z, nrow)))
  #   rownames(z) <- NULL
  # }else{
  #   z$site <- ""
  #   z$site <- factor(z$site)
  # }
  
  z <- fortify(x, type=type)
  
  
  if("y.lwr" %in% names(z) == FALSE & se) {
    warning("invalid se setting, the iNEXT object do not consist confidence interval")
    se <- FALSE
  }else if("y.lwr" %in% names(z) & se) {
    se <- TRUE
  }else{
    se <- FALSE
  }
  
  if(type==1L) {
    #z$x <- z[,1]
    #z$y <- z$qD
    if(!is.null(xlab)) xlab <- ifelse(names(x$DataInfo)[2]=="n", "Number of individuals", "Number of sampling units")
    if(!is.null(ylab)) ylab <- "Species diversity"
    # if(se){
    #   z$y.lwr <- z[,5]
    #   z$y.upr <- z[,6]
    # }
  }else if(type==2L){
    if(length(unique(z$Order.q))>1){
      # z <- subset(z, Order.q==unique(z$Order.q)[1])
      z <- z[z$Order.q==unique(z$Order.q)[1],]
    }
    # z$x <- z[,1]
    # z$y <- z$SC
    if(!is.null(xlab)) xlab <- ifelse(names(x$DataInfo)[2]=="n", "Number of individuals", "Number of sampling units")
    if(!is.null(ylab)) ylab <- "Sample coverage"
    # if(se){
    #   z$y.lwr <- z[,8]
    #   z$y.upr <- z[,9]
    # }
  }else if(type==3L){
    # z$x <- z$SC
    # z$y <- z$qD
    if(!is.null(xlab)) xlab <- "Sample coverage"
    if(!is.null(ylab)) ylab <- "Species diversity"
    # if(se){
    #   z$y.lwr <- z[,5]
    #   z$y.upr <- z[,6]
    # }
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
  
  SITE <- unique(z$Assemblage)
  ORDER <- unique(z$Order.q)
  
  if(is.null(col)){
    col <- gg_color_hue(length(SITE))
  }else{
    col <- rep(col,length(SITE))[1:length(SITE)]
  }
  pch <- (16+1:length(SITE))%%25
  
  for(j in 1:length(ORDER)){
    if(se==TRUE){
      # tmp.sub <- subset(z, Order.q==ORDER[j])
      tmp.sub <- z[z$Order.q==ORDER[j],]
      tmp.j <- data.frame(Assemblage=tmp.sub$Assemblage, Order.q=tmp.sub$Order.q,
                          Method=tmp.sub$Method, 
                          x=tmp.sub$x, y=tmp.sub$y,
                          y.lwr=tmp.sub$y.lwr, y.upr=tmp.sub$y.upr)
      
      plot(y.upr~x, data=tmp.j, type="n", xlab="", ylab="", ...)
    }else{
      # tmp.sub <- subset(z, Order.q==ORDER[j])
      tmp.sub <- z[z$Order.q==ORDER[j],]
      
      tmp.j <- data.frame(Assemblage=tmp.sub$Assemblage, Order.q=tmp.sub$Order.q,
                          Method=tmp.sub$Method, 
                          x=tmp.sub$x, y=tmp.sub$y)
      
      plot(y~x, data=tmp.j, type="n", xlab="", ylab="", ...)
    }
    
    for(i in 1:length(SITE)){
      # tmp <- subset(tmp.j, Assemblage==SITE[i])
      tmp <- tmp.j[tmp.j$Assemblage==SITE[i],]
      if(se==TRUE){
        conf.reg(x=tmp$x, LCL=tmp$y.lwr, UCL=tmp$y.upr, border=NA, col=adjustcolor(col[i], 0.25))
      }
      # lines(y~x, data=subset(tmp, Method=="Rarefaction"), lty=1, lwd=2, col=col[i])
      lines(y~x, data=tmp[tmp$Method=="Rarefaction",], lty=1, lwd=2, col=col[i])
      # lines(y~x, data=subset(tmp, Method=="Extrapolation"), lty=2, lwd=2, col=col[i])
      lines(y~x, data=tmp[tmp$Method=="Extrapolation",], lty=2, lwd=2, col=col[i])
      # points(y~x, data=subset(tmp, Method=="Observed"), pch=pch[i], cex=2, col=col[i])
      points(y~x, data=tmp[tmp$Method=="Observed",], pch=pch[i], cex=2, col=col[i])
      
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
  order.n <- paste(unique(x$iNextEst$size_based$Order.q), collapse = ", ")
  cat("Compare ", site.n, " assemblages with Hill number order q = ", order.n,".\n", sep="")
  cat("$class: iNEXT\n\n")
  cat("$DataInfo: basic data information\n")
  print(x$DataInfo)
  cat("\n")
  cat("$iNextEst: diversity estimates with rarefied and extrapolated samples.\n")
  cat("$size_based (LCL and UCL are obtained for fixed size.)\n")
  cat("\n")
  if(inherits(x$iNextEst, "data.frame")){
    y <- x$iNextEst
    m <- quantile(y[,1], type = 1)
    res <- y[y[,1]%in%m,]
  }else{
    res <- lapply((x$iNextEst), function(y){
      Assemblages <- unique(x$iNextEst$size_based$Assemblage)
      tmp <- lapply(1:length(Assemblages),function(i){
        # y_each <- subset(y, Assemblage==Assemblages[i])
        y_each <- y[y$Assemblage==Assemblages[i],]
        m <- quantile(y_each[,2], type = 1)
        y_each[y_each[,2]%in%m,]
      })
      do.call(rbind,tmp)
    })
  }
  print(res[[1]])
  cat("\n")
  cat("NOTE: The above output only shows five estimates for each assemblage; call iNEXT.object$iNextEst$size_based to view complete output.\n")
  cat("\n")
  cat("$coverage_based (LCL and UCL are obtained for fixed coverage; interval length is wider due to varying size in bootstraps.)\n")
  cat("\n")
  print(res[[2]])
  cat("\n")
  cat("NOTE: The above output only shows five estimates for each assemblage; call iNEXT.object$iNextEst$coverage_based to view complete output.\n")
  cat("\n")
  cat("$AsyEst: asymptotic diversity estimates along with related statistics.\n")
  print(x$AsyEst)
  return(invisible())
}
