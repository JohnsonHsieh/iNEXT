#' iNterpolation and EXTrapolation of Hill number
#' 
#' \code{iNEXT_}: Interpolation and extrapolation of Hill number with order q
#' 
#' @param x a matrix, data.frame (species by sites), or list of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units in each column or list. 
#' @param q a numerical vector of the order of Hill number.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
# @param rowsum a logical variable to check if the input object is raw data (species by sites matrix, \code{rowsum=FALSE}) or iNEXT default input (abundance counts or incidence frequencies, \code{rowsum=TRUE}).
#' @param size an integer vector of sample sizes (number of individuals or sampling units) for which diversity estimates will be computed. 
#' If NULL, then diversity estimates will be computed for those sample sizes determined by the specified/default \code{endpoint} and \code{knots} .
#' @param endpoint an integer specifying the sample size that is the \code{endpoint} for rarefaction/extrapolation. 
#' If NULL, then \code{endpoint} \code{=} double reference sample size.
#' @param knots an integer specifying the number of equally-spaced \code{knots} (say K, default is 40) between size 1 and the \code{endpoint};
#' each knot represents a particular sample size for which diversity estimate will be calculated.  
#' If the \code{endpoint} is smaller than the reference sample size, then \code{iNEXT()} computes only the rarefaction esimates for approximately K evenly spaced \code{knots}. 
#' If the \code{endpoint} is larger than the reference sample size, then \code{iNEXT()} computes rarefaction estimates for approximately K/2 evenly spaced \code{knots} between sample size 1 and the reference sample size, and computes extrapolation estimates for approximately K/2 evenly spaced \code{knots} between the reference sample size and the \code{endpoint}.
#' @param se a logical variable to calculate the bootstrap standard error and \code{conf} confidence interval.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @param nboot an integer specifying the number of replications.
#' @importFrom reshape2 dcast
#' @return a list of three objects: \code{$DataInfo} for summarizing data information; 
#' \code{$iNextEst} for showing diversity estimates for rarefied and extrapolated samples along with related statistics;
#' and \code{$AsyEst} for showing asymptotic diversity estimates along with related statistics.  
#' @examples
#' ## example for abundance based data (list of vector)
#' data(spider)
#' out1 <- iNEXT_(spider, q=0, datatype="abundance")
#' out1$DataInfo # showing basic data information.
#' out1$AsyEst # showing asymptotic diversity estimates.
#' out1$iNextEst # showing diversity estimates with rarefied and extrapolated.
#' 
#' ## example for abundance based data (data.frame)
#' data(bird)
#' out2 <- iNEXT_(bird, q=0, datatype="abundance")
#' ggiNEXT_(out2)
#' 
#' \dontrun{
#' ## example for incidence frequencies based data (list of data.frame)
#' data(ant)
#' t <- round(seq(10, 500, length.out=20))
#' out3 <- iNEXT_(ant$h500m, q=1, datatype="incidence_freq", size=t, se=FALSE)
#' out3$iNextEst
#' }
#' @export
#' 
iNEXT_ <- function(x, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)
{
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  class_x <- class(x)[1]
  
  if(datatype == "incidence"){
    stop('datatype="incidence" was no longer supported after v2.0.8, 
         please try datatype="incidence_freq".')  
  }
  if(datatype=="incidence_freq") datatype <- "incidence"
  
  if(datatype=="incidence_freq") datatype <- "incidence"
  
  if(datatype=="incidence_raw"){
    if(class_x=="list"){
      x <- lapply(x, as.incfreq)
    }else{
      x <- as.incfreq(x)
    }
    datatype <- "incidence"
  }
  #   if(rowsum==FALSE){
  #     if(datatype=="abundance"){
  #       if(class_x =="list"){
  #         x <- lapply(x, as.abucount)
  #       }else{
  #         x <- as.abucount(x)
  #       }
  #     }else if(datatype=="incidence"){
  #       if(class_x =="list"){
  #         x <- lapply(x, as.incfreq)
  #       }else{
  #         x <- as.incfreq(x)
  #       }
  #     }
  #   }
  
  Fun <- function(x, q){
    x <- as.numeric(unlist(x))
    if(datatype == "abundance"){
      if(sum(x)==0) stop("Zero abundance counts in one or more sample sites")
      out <- iNEXT.Ind_(Spec=x, q=q, m=size, endpoint=ifelse(is.null(endpoint), 2*sum(x), endpoint), knots=knots, se=se, nboot=nboot, conf=conf)
    }
    if(datatype == "incidence"){
      t <- x[1]
      y <- x[-1]
      if(t>sum(y)){
        warning("Insufficient data to provide reliable estimators and associated s.e.") 
      }
      if(sum(x)==0) stop("Zero incidence frequencies in one or more sample sites")
      
      out <- iNEXT.Sam_(Spec=x, q=q, t=size, endpoint=ifelse(is.null(endpoint), 2*max(x), endpoint), knots=knots, se=se, nboot=nboot, conf=conf)  
    }
    out
  }
  
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a postive value/vector of numeric object")
  if(min(q) < 0){
    warning("ambigous of order q, we only compute postive q")
    q <- q[q >= 0]
  }
  
  z <- qnorm(1-(1-0.95)/2)
  if(class_x=="numeric" | class_x=="integer" | class_x=="double"){
    out <- Fun(x,q)
    index <- AsymDiv(x = x,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq')
                     ,nboot = 100,conf = 0.95,method = 'Both')
    LCL <- index$qD.LCL[index$method=='Estimated']
    UCL <- index$qD.UCL[index$method=='Estimated']
    index <- dcast(index,formula = order~method,value.var = 'qD')
    index <- cbind(index[,-1],se = (UCL - index$Estimated)/z,LCL,UCL)
    index$LCL[index$LCL<index$Empirical & index$order==0] <- index$Empirical[index$LCL<index$Empirical & index$order==0]
    colnames(index) <- c("Observed","Estimator","Est_s.e.","95% Lower","95% Upper")
    # index <- rbind(as.matrix(ChaoSpecies(x, datatype, conf)), 
    #                as.matrix(ChaoEntropy(x, datatype, transform=TRUE, conf)),
    #                as.matrix(EstSimpson(x, datatype, transform=TRUE, conf)))
    rownames(index) <- c("Species Richness", "Shannon diversity", "Simpson diversity")
    
    
  }else if(class_x=="matrix" | class_x=="data.frame"){
    out <- apply(as.matrix(x), 2, function(x){
      tmp <- Fun(x,q)
      tmp
    })
    index <- AsymDiv(x = x,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq'),nboot = 100,conf = 0.95,method = 'Both')
    LCL <- index$qD.LCL[index$method=='Estimated']
    UCL <- index$qD.UCL[index$method=='Estimated']
    index <- dcast(index,formula = Site+order~method,value.var = 'qD')
    index <- cbind(index,se = (UCL - index$Estimated)/z,LCL,UCL)
    index$LCL[index$LCL<index$Empirical & index$order==0] <- index$Empirical[index$LCL<index$Empirical & index$order==0]
    index$order <- c('Species richness','Shannon diversity','Simpson diversity')
    # arr <- array(0, dim = c(3, 5, ncol(x)))
    # arr[1,,] <- t(as.matrix(ChaoSpecies(x, datatype, conf)))
    # arr[2,,] <- t(as.matrix(ChaoEntropy(x, datatype, transform=TRUE, conf)))
    # arr[3,,] <- t(as.matrix(EstSimpson(x, datatype, transform=TRUE, conf)))  
    # dimnames(arr)[[3]] <- names(x)
    # dimnames(arr)[[1]] <- c("Species richness", "Shannon diversity", "Simpson diversity")
    # dimnames(arr)[[2]] <- c("Observed", "Estimator", "Est_s.e.", "Lower_CI", "Upper_CI")
    # index <- ftable(arr, row.vars = c(3,1))
    # index <- dcast(as.data.frame(index), formula = Var1+Var2~Var3, value.var = "Freq")
    colnames(index) <- c("Site", "Diversity", "Observed", "Estimator", "s.e.", "LCL", "UCL")
    
    
  }else if(class_x=="list"){
    out <- lapply(x, function(x) {
      tmp <- Fun(x,q)
      tmp
    })
    index <- AsymDiv(x = x,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq'),nboot = 100,conf = 0.95,method = 'Both')
    LCL <- index$qD.LCL[index$method=='Estimated']
    UCL <- index$qD.UCL[index$method=='Estimated']
    index <- dcast(index,formula = Site+order~method,value.var = 'qD')
    index <- cbind(index,se = (UCL - index$Estimated)/z,LCL,UCL)
    index$LCL[index$LCL<index$Empirical & index$order==0] <- index$Empirical[index$LCL<index$Empirical & index$order==0]
    index$order <- c('Species richness','Shannon diversity','Simpson diversity')
    # arr <- array(0, dim = c(3, 5, length(x)))
    # arr[1,,] <- t(as.matrix(ChaoSpecies(x, datatype, conf)))
    # arr[2,,] <- t(as.matrix(ChaoEntropy(x, datatype, transform=TRUE, conf)))
    # arr[3,,] <- t(as.matrix(EstSimpson(x, datatype, transform=TRUE, conf)))  
    # dimnames(arr)[[3]] <- names(x)
    # dimnames(arr)[[1]] <- c("Species richness", "Shannon diversity", "Simpson diversity")
    # dimnames(arr)[[2]] <- c("Observed", "Estimator", "Est_s.e.", "Lower_CI", "Upper_CI")
    # index <- ftable(arr, row.vars = c(3,1))
    # index <- dcast(as.data.frame(index), formula = Var1+Var2~Var3, value.var = "Freq")
    colnames(index) <- c("Site", "Diversity", "Observed", "Estimator", "s.e.", "LCL", "UCL")
  }else{
    stop("invalid class of x, x should be a object of numeric, matrix, data.frame, or list")
  }
  
  info <- DataInfo(x, datatype)
  
  
  z <- list("DataInfo"=info, "iNextEst"=out, "AsyEst"=index)
  class(z) <- c("iNEXT_")
  return(z)
}

# iNterpolation and EXTrapolation of abundance-based Hill number
# 
# \code{iNEXT.Ind_} Estimation of interpolation and extrapolation of abundance-based Hill number with order q
# 
# @param Spec a vector of species abundances
# @param q a numerical vector of the order of Hill number
# @param m a integer vector of rarefaction/extrapolation sample size, default is NULL. If m is not be specified, then the program will compute sample units due to endpoint and knots.
# @param endpoint a integer of sample size that is the endpoint for rarefaction/extrapolation. Default is double the original sample size.
# @param knots a number of knots of computation, default is 40
# @param se calculate bootstrap standard error and 95% confidence interval; default is TRUE
# @param nboot the number of bootstrap resampling times, default is 200
# @return a list of interpolation and extrapolation Hill number with specific order q (qD) and sample coverage (SC)
# @seealso \code{\link{iNEXT.Sam_}}
# @examples
# data(spider)
# # q = 0 with specific endpoint
# iNEXT.Ind_(spider$Girdled, q=0, endpoint=500)
# # q = 1 with specific sample size m and don't calculate standard error
# iNEXT.Ind_(spider$Girdled, q=1, m=c(1, 10, 20, 50, 100, 200, 400, 600), se=FALSE)
iNEXT.Ind_ <- function(Spec, q=0, m=NULL, endpoint=2*sum(Spec), knots=40, se=TRUE, nboot=200, conf=0.95)
{
  
  n <- sum(Spec)		  	#sample size
  if(is.null(m)) {
    if(endpoint <= n) {
      m <- floor(seq(1, endpoint, length=floor(knots)))
    } else {
      m <- c(floor(seq(1, sum(Spec)-1, length.out=floor(knots/2)-1)), sum(Spec), floor(seq(sum(Spec)+1, to=endpoint, length.out=floor(knots/2))))
    }
    m <- c(1, m[-1])
  } else if(is.null(m)==FALSE) {	
    if(max(m)>n & length(m[m==n])==0)  m <- c(m, n-1, n, n+1)
    m <- sort(m)
  }
  
  Dq.hat <- TD.m.est(Spec,m,q)
  C.hat <- Chat.Ind(Spec, m)
  
  if(se==TRUE & nboot > 0 & length(Spec) > 1) {
    Prob.hat <- EstiBootComm.Ind(Spec)
    Abun.Mat <- rmultinom(nboot, n, Prob.hat)
    
    error <-  qnorm(1-(1-conf)/2) * apply(apply(Abun.Mat, 2, function(x) TD.m.est(x, m, q)), 1, sd, na.rm=TRUE)
    left  <- Dq.hat - error
    right <- Dq.hat + error
    
    error.C <-  qnorm(1-(1-conf)/2) * apply(apply(Abun.Mat, 2, function(x) Chat.Ind(x, m)), 1, sd, na.rm=TRUE)
    left.C  <- C.hat - error.C
    right.C <- C.hat + error.C
    out <- cbind("m"=rep(m,length(q)), "qD"=Dq.hat, "qD.LCL"=left, "qD.UCL"=right, 
                 "SC"=rep(C.hat,length(q)), "SC.LCL"=left.C, "SC.UCL"=right.C)
  } else {
    out <- cbind("m"=rep(m,length(q)), "qD"=Dq.hat, "SC"=rep(C.hat,length(q)))
  }
  out <- data.frame(out)
  out$method <- ifelse(out$m<n, "interpolated", ifelse(out$m==n, "observed", "extrapolated"))
  out$order <- rep(q,each = length(m))
  id <- match(c("m", "method", "order", "qD", "qD.LCL", "qD.UCL", "SC", "SC.LCL", "SC.UCL"), names(out), nomatch = 0)
  out <- out[, id]
  return(out)
}

###############################################
# iNterpolation and EXTrapolation of incidence-based Hill number
# 
# \code{iNEXT.Sam_} Estimation of interpolation and extrapolation of incidence-based Hill number with order q
# 
# @param Spec a vector of species incidence-based frequency, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
# @param q a numerical vector of the order of Hill number
# @param t a integer vector of rarefaction/extrapolation sample size, default is NULL. If m is not be specified, then the program will compute sample units due to endpoint and knots.
# @param endpoint a integer of sample size that is the endpoint for rarefaction/extrapolation. Default is double the original sample size.
# @param knots a number of knots of computation, default is 40
# @param se calculate bootstrap standard error and 95% confidence interval; default is TRUE
# @param nboot the number of bootstrap resampling times, default is 200
# @return a list of interpolation and extrapolation Hill number with specific order q (qD) and sample coverage (SC)
# @seealso \code{\link{iNEXT.Ind_}}
# @examples
# data(ant)
# # q = 0 with specific endpoint
# iNEXT.Sam_(ant$h50m, q=0, endpoint=100)
# # q = 1 with specific sample size m and don't calculate standard error
# iNEXT.Sam_(ant$h500m, q=1, t=round(seq(10, 500, length.out=20)), se=FALSE)
iNEXT.Sam_ <- function(Spec, t=NULL, q=0, endpoint=2*max(Spec), knots=40, se=TRUE, nboot=200, conf=0.95)
{
  
  if(which.max(Spec)!=1) 
    stop("invalid data structure!, first element should be number of sampling units")
  
  nT <- Spec[1]
  if(is.null(t)) {
    if(endpoint <= nT) {
      t <- floor(seq(1, endpoint, length.out=floor(knots)))
    } else {
      t <- c(floor(seq(1, nT-1, length.out=floor(knots/2)-1)), nT, floor(seq(nT+1, to=endpoint, length.out=floor(knots/2))))
    }
    t <- c(1, t[-1])
  } else if(is.null(t)==FALSE) {	
    if(max(t)>nT & length(t[t==nT])==0)  t <- c(t, nT-1, nT, nT+1)
    t <- sort(t)
  }
  
  Dq.hat <- TD.m.est_inc(Spec,t,q)
  C.hat <- Chat.Sam(Spec, t)
  
  if(se==TRUE & nboot > 0 & length(Spec) > 2){
    Prob.hat <- EstiBootComm.Sam(Spec)
    Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
    Abun.Mat <- matrix(c(rbind(nT, Abun.Mat)),ncol=nboot)
    tmp <- which(colSums(Abun.Mat)==nT)
    if(length(tmp)>0) Abun.Mat <- Abun.Mat[,-tmp]
    if(ncol(Abun.Mat)==0){
      out <- cbind("t"=t, "qD"=Dq.hat, "SC"=C.hat)
      warning("Insufficient data to compute bootstrap s.e.")
    }else{		
      error <-  qnorm(1-(1-conf)/2) * apply(apply(Abun.Mat, 2, function(y) TD.m.est_inc(y,t,q)), 1, sd, na.rm=TRUE)
      left  <- Dq.hat - error
      right <- Dq.hat + error
      left[left<=0] <- 0
      
      error.C <-  qnorm(1-(1-conf)/2) * apply(apply(Abun.Mat, 2, function(y) Chat.Sam(y, t)), 1, sd, na.rm=TRUE)
      left.C  <- C.hat - error.C
      right.C <- C.hat + error.C
      left.C[left.C<=0] <- 0
      right.C[right.C>=1] <- 1
      
      out <- cbind("t"=rep(t,length(q)), "qD"=Dq.hat, "qD.LCL"=left, "qD.UCL"=right, 
                   "SC"=rep(C.hat,length(q)), "SC.LCL"=left.C, "SC.UCL"=right.C)
    }
  }else {
    out <- cbind("t"=rep(t,length(q)), "qD"=Dq.hat, "SC"=rep(C.hat,length(q)))
  }
  out <- data.frame(out)
  out$method <- ifelse(out$t<nT, "interpolated", ifelse(out$t==nT, "observed", "extrapolated"))
  out$order <- rep(q,each = length(t))
  id <- match(c("t", "method", "order", "qD", "qD.LCL", "qD.UCL", "SC", "SC.LCL", "SC.UCL"), names(out), nomatch = 0)
  out <- out[, id]
  return(out)
}

#=====For plotting outcome=====
#
#
###############################################
#' ggplot2 extension for an iNEXT_ object
#' 
#' \code{ggiNEXT_}: the \code{\link[ggplot2]{ggplot}} extension for \code{\link{iNEXT_}} Object to plot sample-size- and coverage-based rarefaction/extrapolation curves along with a bridging sample completeness curve
#' @param x an \code{iNEXT_} object computed by \code{\link{iNEXT_}}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).                 
#' @param se a logical variable to display confidence interval around the estimated sampling curve.
#' @param facet.var create a separate plot for each value of a specified variable: 
#'  no separation \cr (\code{facet.var="none"}); 
#'  a separate plot for each diversity order (\code{facet.var="order"}); 
#'  a separate plot for each site (\code{facet.var="site"}); 
#'  a separate plot for each combination of order x site (\code{facet.var="both"}).              
#' @param color.var create curves in different colors for values of a specified variable:
#'  all curves are in the same color (\code{color.var="none"}); 
#'  use different colors for diversity orders (\code{color.var="order"}); 
#'  use different colors for sites (\code{color.var="site"}); 
#'  use different colors for combinations of order x site (\code{color.var="both"}).  
#' @param grey a logical variable to display grey and white ggplot2 theme. 
#' @param ... other arguments passed on to methods. Not currently used.
#' @return a ggplot2 object
#' @examples
#' data(spider)
#' # single-assemblage abundance data
#' out1 <- iNEXT_(spider$Girdled, q=0, datatype="abundance")
#' ggiNEXT_(x=out1, type=1)
#' ggiNEXT_(x=out1, type=2)
#' ggiNEXT_(x=out1, type=3)
#' 
#'\dontrun{
#' # single-assemblage incidence data with three orders q
#' data(ant)
#' size <- round(seq(10, 500, length.out=20))
#' y <- iNEXT_(ant$h500m, q=c(0,1,2), datatype="incidence_freq", size=size, se=FALSE)
#' ggiNEXT_(y, se=FALSE, color.var="order")
#' 
#' # multiple-assemblage abundance data with three orders q
#' z <- iNEXT_(spider, q=c(0,1,2), datatype="abundance")
#' ggiNEXT_(z, facet.var="site", color.var="order")
#' ggiNEXT_(z, facet.var="both", color.var="both")
#'}
#' @export
#' 
ggiNEXT_ <- function(x, type=1, se=TRUE, facet.var="none", color.var="site", grey=FALSE){  
  UseMethod("ggiNEXT_", x)
}

#' @export
#' @rdname ggiNEXT_
ggiNEXT_.iNEXT_ <- function(x, type=1, se=TRUE, facet.var="none", color.var="site", grey=FALSE){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#330066", "#CC79A7",  "#0072B2", "#D55E00"))
  TYPE <-  c(1, 2, 3)
  SPLIT <- c("none", "order", "site", "both")
  if(is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1)
    stop("invalid plot type")
  if(is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var, SPLIT) == -1)
    stop("invalid facet variable")
  if(is.na(pmatch(color.var, SPLIT)) | pmatch(color.var, SPLIT) == -1)
    stop("invalid color variable")
  
  type <- pmatch(type, 1:3)
  facet.var <- match.arg(facet.var, SPLIT)
  color.var <- match.arg(color.var, SPLIT)
  if(facet.var=="order") color.var <- "site"
  if(facet.var=="site") color.var <- "order"
  
  options(warn = -1)
  z <- fortify(x, type=type)
  options(warn = 0)
  if(ncol(z) ==7) {se <- FALSE}
  datatype <- unique(z$datatype)
  if(color.var=="none"){
    if(levels(factor(z$order))>1 & "site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT_ object consists multiple sites and orders, change setting as both")
      color.var <- "both"
      z$col <- z$shape <- paste(z$site, z$order, sep="-")
      
    }else if("site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT_ object consists multiple orders, change setting as order")
      color.var <- "site"
      z$col <- z$shape <- z$site
    }else if(levels(factor(z$order))>1){
      warning("invalid color.var setting, the iNEXT_ object consists multiple sites, change setting as site")
      color.var <- "order"
      z$col <- z$shape <- factor(z$order)
    }else{
      z$col <- z$shape <- rep(1, nrow(z))
    }
  }else if(color.var=="order"){     
    z$col <- z$shape <- factor(z$order)
  }else if(color.var=="site"){
    if(!"site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT_ object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- z$site
  }else if(color.var=="both"){
    if(!"site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT_ object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- paste(z$site, z$order, sep="-")
  }
  zz=z
  z$method[z$method=="observed"]="interpolated"
  z$lty <- z$lty <- factor(z$method, levels=unique(c("interpolated", "extrapolated"),
                                                   c("interpolation", "interpolation", "extrapolation")))
  z$col <- factor(z$col)
  data.sub <- zz[which(zz$method=="observed"),]
  
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
  }
  else if(type==3L) {
    g <- g + labs(x="Sample coverage", y="Species diversity")
  }
  else {
    g <- g + labs(x="Number of sampling units", y="Species diversity")
    if(datatype=="abundance") g <- g + labs(x="Number of individuals", y="Species diversity")
  }
  
  if(se)
    g <- g + geom_ribbon(aes_string(ymin="y.lwr", ymax="y.upr", fill="factor(col)", colour="NULL"), alpha=0.2)+
    scale_fill_manual(values=cbPalette)
  
  
  if(facet.var=="order"){
    if(length(levels(factor(z$order))) == 1 & type!=2){
      warning("invalid facet.var setting, the iNEXT_ object do not consist multiple orders.")      
    }else{
      odr_grp <- as_labeller(c(`0` = "q = 0", `1` = "q = 1",`2` = "q = 2")) 
      g <- g + facet_wrap(~order, nrow=1, labeller = odr_grp)
      if(color.var=="both"){
        g <- g + guides(colour=guide_legend(title="Guides", ncol=length(levels(factor(z$order))), byrow=TRUE),
                        fill=guide_legend(title="Guides"))
      }
      if(type==2){
        g <- g + theme(strip.background = element_blank(),strip.text.x = element_blank())
        
      }
    }
  }
  
  if(facet.var=="site"){
    if(!"site"%in%names(z)) {
      warning("invalid facet.var setting, the iNEXT_ object do not consist multiple sites.")
    }else{
      g <- g + facet_wrap(~site, nrow=1)
      if(color.var=="both"){
        g <- g + guides(colour=guide_legend(title="Guides", nrow=length(levels(factor(z$order)))),
                        fill=guide_legend(title="Guides"))
      }
    }
  }
  
  if(facet.var=="both"){
    if(length(levels(factor(z$order))) == 1 | !"site"%in%names(z)){
      warning("invalid facet.var setting, the iNEXT_ object do not consist multiple sites or orders.")
    }else{
      odr_grp <- as_labeller(c(`0` = "q = 0", `1` = "q = 1",`2` = "q = 2")) 
      g <- g + facet_wrap(site~order,labeller = labeller(order = odr_grp)) 
      if(color.var=="both"){
        g <- g +  guides(colour=guide_legend(title="Guides", nrow=length(levels(factor(z$site))), byrow=TRUE),
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
#' @rdname ggiNEXT_
ggiNEXT_.default <- function(x, ...){
  stop(
    "iNEXT_ doesn't know how to deal with data of class ",
    paste(class(x), collapse = "/"),
    call. = FALSE
  )
}

#' Fortify method for classes from the iNEXT package.
#'
#' @name fortify.iNEXT_
#' @param model \code{iNEXT_} to convert into a dataframe.
#' @param data not used by this method
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).                 
#' @param ... not used by this method
#' @import ggplot2
#' @export
#' @examples
#' data(spider)
#' # single-assemblage abundance data
#' out1 <- iNEXT_(spider$Girdled, q=0, datatype="abundance")
#' ggplot2::fortify(out1, type=1)
fortify.iNEXT_ <- function(model, data = model$iNextEst, type = 1, ...) {
  datatype <- ifelse(names(model$DataInfo)[2]=="n","abundance","incidence")
  z <- data
  if(class(z) == "list"){
    z <- data.frame(do.call("rbind", z), site=rep(names(z), sapply(z, nrow)))
    rownames(z) <- NULL
  }else{
    z$site <- ""
  }
  
  if(ncol(z)==6) {
    warning("invalid se setting, the iNEXT_ object do not consist confidence interval")
    se <- FALSE
  }else if(ncol(z)>6) {
    se <- TRUE
  }
  
  if(type==1L) {
    z$x <- z[,1]
    z$y <- z$qD
    if(se){
      z$y.lwr <- z[,5]
      z$y.upr <- z[,6]
    }
  }else if(type==2L){
    if(length(unique(z$order))>1){
      z <- subset(z, order==unique(z$order)[1])
    }
    z$x <- z[,1]
    z$y <- z$SC
    if(se){
      z$y.lwr <- z[,8]
      z$y.upr <- z[,9]
    }
  }else if(type==3L){
    z$x <- z$SC
    z$y <- z$qD
    if(se){
      z$y.lwr <- z[,5]
      z$y.upr <- z[,6]
    }
  }
  z$datatype <- datatype
  z$plottype <- type
  if(se){
    data <- z[,c("datatype","plottype","site","method","order","x","y","y.lwr","y.upr")]
  }else{
    data <- z[,c("datatype","plottype","site","method","order","x","y")]
  }
  data
}

###############################################
#' Plotting iNEXT_ object
#' 
#' \code{plot.iNEXT_}: Plotting method for objects inheriting from class "iNEXT_"
#' @param x an \code{iNEXT_} object computed by \code{\link{iNEXT_}}.
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
#' out1 <- iNEXT_(spider$Girdled, q=0, datatype="abundance")
#' plot(x=out1, type=1)
#' plot(x=out1, type=2)
#' plot(x=out1, type=3)
#' 

#' @export
plot.iNEXT_ <- function(x, type=1, se=TRUE, show.legend=TRUE, show.main=TRUE, col=NULL,...){
  
  if(class(x) != "iNEXT_")
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
      z$y.lwr <- z[,5]
      z$y.upr <- z[,6]
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
      z$y.lwr <- z[,8]
      z$y.upr <- z[,9]
    }
  }else if(type==3L){
    z$x <- z$SC
    z$y <- z$qD
    if(!is.null(xlab)) xlab <- "Sample coverage"
    if(!is.null(ylab)) ylab <- "Species diversity"
    if(se){
      z$y.lwr <- z[,5]
      z$y.upr <- z[,6]
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


#' Printing iNEXT_ object
#' 
#' \code{print.iNEXT_}: Print method for objects inheriting from class "iNEXT_"
#' @param x an \code{iNEXT_} object computed by \code{\link{iNEXT_}}.
#' @param ... additional arguments.
#' @export
print.iNEXT_ <- function(x, ...){
  site.n <- nrow(x$DataInfo)
  order.n <- ifelse(site.n > 1, 
                    paste(unique(x$iNextEst[[1]]$order), collapse = ", "),
                    paste(unique(x$iNextEst$order), collapse = ", "))
  cat("Compare ", site.n, " assemblages with Hill number order q = ", order.n,".\n", sep="")
  cat("$class: iNEXT_\n\n")
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
  cat("NOTE: Only show five estimates, call iNEXT.object$iNextEst. to show complete output.\n")
  return(invisible())
}

