invChat.Ind <- function(x, C, conf=NULL){
  m <- NULL # no visible binding for global variable 'm'
  n <- sum(x)
  refC <- Chat.Ind(x,n)
  f <- function(m, C) abs(Chat.Ind(x,m)-C)
  if(refC > C)
  {
     opt <- optimize(f, C=C, lower=0, upper=sum(x))
     mm <- opt$minimum
     mm <- round(mm)
  }
  if(refC <= C)
  {
     f1 <- sum(x==1)
     f2 <- sum(x==2)
     if(f1>0 & f2>0){A <- (n-1)*f1/((n-1)*f1+2*f2)}
     if(f1>1 & f2==0){A <- (n-1)*(f1-1)/((n-1)*(f1-1)+2)}
     if(f1==1 & f2==0){A <- 1}
     if(f1==0 & f2==0){A <- 1}
     mm <- (log(n/f1)+log(1-C))/log(A)-1
     mm <- n+mm
     mm <- round(mm)
  }
  if(mm > 2*n) 
	warning("The maximum size of the extrapolation exceeds double reference sample size, the results for q = 0 may be subject to large prediction bias.")
  method <- ifelse(mm<n, "interpolated", ifelse(mm==n, "observed", "extrapolated"))
  # out <- data.frame(m=mm, method=method, 
  #                   "SC"=round(Chat.Ind(x,mm),4), 
  #                   "q = 0"=round(Dqhat.Ind(x,0,mm),3),
  #                   "q = 1"=round(Dqhat.Ind(x,1,mm),3),
  #                   "q = 2"=round(Dqhat.Ind(x,2,mm),3))
  # colnames(out) <- c("m", "method", "SC", "q = 0", "q = 1", "q = 2")
  # out
  
  if(is.null(conf)){
    out <- data.frame(m=mm, method=method, 
                      SamCov=round(Chat.Ind(x,mm),3),
                      SpeRic=round(Dqhat.Ind(x,0,mm),3),
                      ShaDiv=round(Dqhat.Ind(x,1,mm),3),
                      SimDiv=round(Dqhat.Ind(x,2,mm),3))
    colnames(out) <- c("m", "method", "SC", "q = 0", "q = 1", "q = 2")
  }else{
    tmp0 <- iNEXT.Ind(x, q=0, m=c(1,mm), se = TRUE, conf = conf)
    tmp1 <- iNEXT.Ind(x, q=1, m=c(1,mm), se = TRUE, conf = conf)
    tmp2 <- iNEXT.Ind(x, q=2, m=c(1,mm), se = TRUE, conf = conf)
    tmp <- subset(rbind(tmp0,tmp1,tmp2), m==mm)
    out <- tmp[,c(1,2,3,7,4,5,6)]
    out[,4:7] <- round(out[,4:7],3)
  }
  out
}

invChat.Sam <- function(x, C, conf=NULL){
  m <- NULL # no visible binding for global variable 'm'
  n <- max(x)
  refC <- Chat.Sam(x,n)
  f <- function(m, C) abs(Chat.Sam(x,m)-C)
  if(refC > C)
  {
     opt <- optimize(f, C=C, lower=0, upper=max(x))
     mm <- opt$minimum
     mm <- round(mm)
  }
  if(refC <= C)
  {
     f1 <- sum(x==1)
     f2 <- sum(x==2)
	 U <- sum(x)-max(x)
     if(f1>0 & f2>0){A <- (n-1)*f1/((n-1)*f1+2*f2)}
     if(f1>1 & f2==0){A <- (n-1)*(f1-1)/((n-1)*(f1-1)+2)}
     if(f1==1 & f2==0){A <- 1}
     if(f1==0 & f2==0){A <- 1}
     mm <- (log(U/f1)+log(1-C))/log(A)-1
     mm <- n+mm
     mm <- round(mm)
  }
  if(mm > 2*n) 
	warning("The maximum size of the extrapolation exceeds double reference sample size, the results for q = 0 may be subject to large prediction bias.")

  
  if(is.null(conf)){
    method <- ifelse(mm<n, "interpolated", ifelse(mm==n, "observed", "extrapolated"))
    out <- data.frame(t=mm, method=method, 
                      SamCov=round(Chat.Sam(x,mm),3),
                      SpeRic=round(Dqhat.Sam(x,0,mm),3),
                      ShaDiv=round(Dqhat.Sam(x,1,mm),3),
                      SimDiv=round(Dqhat.Sam(x,2,mm),3))
    colnames(out) <- c("t", "method", "SC", "q = 0", "q = 1", "q = 2")
  }else{
    tmp0 <- iNEXT.Sam(x, q=0, t=c(1,mm), se = TRUE, conf = conf)
    tmp1 <- iNEXT.Sam(x, q=1, t=c(1,mm), se = TRUE, conf = conf)
    tmp2 <- iNEXT.Sam(x, q=2, t=c(1,mm), se = TRUE, conf = conf)
    tmp <- subset(rbind(tmp0,tmp1,tmp2), t==mm)
    out <- tmp[,c(1,2,3,7,4,5,6)]
    out[,4:7] <- round(out[,4:7],3)
  }
  out
}


invSize.Ind <- function(x, size, conf=NULL){
  m <- NULL # no visible binding for global variable 'm'
  
  if(is.null(size)){
    size <- sum(x)
  }
  method <- ifelse(size<sum(x), "interpolated", ifelse(size==sum(x), "observed", "extrapolated"))
  if(is.null(conf)){
    out <- data.frame(m=size, method=method, 
                      SamCov=round(Chat.Ind(x,size),3),
                      SpeRic=round(Dqhat.Ind(x,0,size),3),
                      ShaDiv=round(Dqhat.Ind(x,1,size),3),
                      SimDiv=round(Dqhat.Ind(x,2,size),3))
    colnames(out) <- c("m", "method", "SC", "q = 0", "q = 1", "q = 2")
  }else{
    tmp0 <- iNEXT.Ind(x, q=0, m=c(1,size), se = TRUE, conf = conf)
    tmp1 <- iNEXT.Ind(x, q=1, m=c(1,size), se = TRUE, conf = conf)
    tmp2 <- iNEXT.Ind(x, q=2, m=c(1,size), se = TRUE, conf = conf)
    tmp <- subset(rbind(tmp0,tmp1,tmp2), m==size)
    out <- tmp[,c(1,2,3,7,4,5,6)]
    out[,4:7] <- round(out[,4:7],3)
  }
  out
}

invSize.Sam <- function(x, size, conf=NULL){
  m <- NULL # no visible binding for global variable 'm'
  
  if(is.null(size)){
    size <- max(x)
  }
  if(is.null(conf)){
    method <- ifelse(size<max(x), "interpolated", ifelse(size==max(x), "observed", "extrapolated"))
    out <- data.frame(t=size, method=method, 
                      SamCov=round(Chat.Sam(x,size),3),
                      SpeRic=round(Dqhat.Sam(x,0,size),3),
                      ShaDiv=round(Dqhat.Sam(x,1,size),3),
                      SimDiv=round(Dqhat.Sam(x,2,size),3))
    colnames(out) <- c("m", "method", "SC", "q = 0", "q = 1", "q = 2")
    out
  }else{
    tmp0 <- iNEXT.Sam(x, q=0, t=c(1,size), se = TRUE, conf = conf)
    tmp1 <- iNEXT.Sam(x, q=1, t=c(1,size), se = TRUE, conf = conf)
    tmp2 <- iNEXT.Sam(x, q=2, t=c(1,size), se = TRUE, conf = conf)
    tmp <- subset(rbind(tmp0,tmp1,tmp2), t==size)
    out <- tmp[,c(1,2,3,7,4,5,6)]
    out[,4:7] <- round(out[,4:7],3)
  }
  out
}
#
#
###############################################
# Compute species diversity with fixed sample coverage
# 
# \code{invChat} compute species diversity with fixed sample coverage
# @param x a \code{data.frame} or \code{list} for species abundance/incidence frequencies.
# @param datatype the data type of input data. That is individual-based abundance data (\code{datatype = "abundance"}) or sample-based incidence data (\code{datatype = "incidence"}).
# @param C a specific sample coverage to compare, which is between 0 to 1. Default is the minimum of double sample size for all sites.
# @return a \code{data.frame} with fixed sample coverage to compare species diversity.
# @examples
# data(spider)
# incChat(spider, "abundance")
# incChat(spider, "abundance", 0.85)
# 
# @export

invChat <- function(x, datatype="abundance", C=NULL, conf=NULL){
  TYPE <- c("abundance", "incidence")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  class_x <- class(x)[1]
  if(datatype=="abundance"){
    if(class_x=="numeric" | class_x=="integer"){
      if(is.null(C)){
        C <- Chat.Ind(x, 2*sum(x))
      }
      out <- invChat.Ind(x, C, conf)
    }
    else if(class_x=="list"){
      if(is.null(C)){
        C <- min(unlist(lapply(x, function(x) Chat.Ind(x, sum(x)))))
      }
      out <- do.call(rbind, lapply(x, function(x) invChat.Ind(x, C, conf)))
    }else if(class_x=="data.frame" | class_x=="matrix"){
      if(is.null(C)){
        C <- min(unlist(apply(x, 2, function(x) Chat.Ind(x, sum(x)))))
      }
      out <- do.call(rbind, apply(x, 2, function(x) invChat.Ind(x, C, conf)))
    }
  }else if(datatype=="incidence"){
    if(class_x=="numeric" | class_x=="integer"){
      if(is.null(C)){
        C <- Chat.Sam(x, 2*max(x))
      }
      out <- invChat.Sam(x, C, conf)
    }
    else if(class_x=="list"){
      if(is.null(C)){
        C <- min(unlist(lapply(x, function(x) Chat.Sam(x, max(x)))))
      }
      out <- do.call(rbind, lapply(x, function(x) invChat.Sam(x, C, conf)))
    }else if(class_x=="data.frame" | class_x=="matrix"){
      if(is.null(C)){
        C <- min(unlist(apply(x, 2, function(x) Chat.Sam(x, max(x)))))
      }
      out <- do.call(rbind, apply(x, 2, function(x) invChat.Sam(x, C, conf)))
    }
  }
  out
}


invSize <- function(x, datatype="abundance", size=NULL, conf=NULL){
  TYPE <- c("abundance", "incidence")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  class_x <- class(x)[1]
  if(datatype=="abundance"){
    if(class_x=="numeric" | class_x=="integer"){
      if(is.null(size)){
        size <- sum(x)
      }
      out <- invSize.Ind(x, size, conf)
    }
    else if(class_x=="list"){
      if(is.null(size)){
        size <- min(unlist(lapply(x, function(x) sum(x))))
      }
      out <- do.call(rbind, lapply(x, function(x) invSize.Ind(x, size, conf)))
    }else if(class_x=="data.frame" | class_x=="matrix"){
      if(is.null(size)){
        size <- min(unlist(apply(x, 2, function(x) sum(x))))
      }
      out <- do.call(rbind, apply(x, 2, function(x) invSize.Ind(x, size, conf)))
    }
  }else if(datatype=="incidence"){
    
    if(class_x=="numeric" | class_x=="integer"){
      if(is.null(size)){
        size <- max(x)
      }
      out <- invSize.Sam(x, size, conf)
    }
    else if(class_x=="list"){
      if(is.null(size)){
        size <- min(unlist(lapply(x, function(x) max(x))))
        }
      out <- do.call(rbind, lapply(x, function(x) invSize.Sam(x, size, conf)))
    }else if(class_x=="data.frame" | class_x=="matrix"){
      if(is.null(size)){
        size <- min(unlist(apply(x, 2, function(x) max(x))))
      }
      out <- do.call(rbind, apply(x, 2, function(x) invSize.Sam(x, size, conf)))
    }
  }
  out
}

# 
# invChat(spider, datatype="abundance")
# invChat(spider, datatype="abundance", C = 0.820)
# invChat(spider, datatype="abundance", C = 0.923)
# invChat(spider, datatype="abundance", C = 0.900)
# 
# invChat(ant, datatype="incidence", 0.95)
# invChat(ant, datatype="incidence", 0.99)
# 

###############################################
#' Compute species diversity with a particular of sample size/coverage 
#' 
#' \code{estimateD}: computes species diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage.
#' @param x a \code{data.frame} or \code{list} of species abundances or incidence frequencies.\cr 
#' If \code{datatype = "incidence"}, then the first entry of the input data must be total number of sampling units, followed 
#' by species incidence frequencies in each column or list.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param base comparison base: sample-size-based (\code{base="size"}) or coverage-based \cr (\code{base="coverage"}).
#' @param level an value specifying a particular sample size or a number (between 0 and 1) specifying a particular value of sample coverage. 
#' If \code{base="size"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample size among all sites. 
#' If \code{base="coverage"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample coverage among all sites. 
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. Remove C.I. by setting conf=NULL.
#' @return a \code{data.frame} of species diversity table including the sample size, sample coverage,
#' method (rarefaction or extrapolation), and diversity estimates with q = 0, 1, and 2 for the user-specified sample size or sample coverage.
#' @examples
#' \dontrun{
#' data(spider)
#' estimateD(spider, "abundance", base="size", level=NULL, conf=0.95)
#' estimateD(spider, "abundance", base="coverage", level=NULL, conf=0.95)
#' }
#' data(ant)
#' estimateD(ant, "incidence_freq", base="coverage", level=0.985, conf=NULL)
#' @export
estimateD <- function(x, datatype="abundance", base="size", level=NULL, conf=0.95){
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  #TYPE <- c("abundance", "incidence")
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
  
  if(datatype=="incidence_raw"){
    if(class_x=="data.frame" | class_x=="matrix") x <- as.incfreq(x)
    else if(class_x=="list") x <- lapply(x, as.incfreq)
    datatype <- "incidence"
  }
    
  BASE <- c("size", "coverage")
  if(is.na(pmatch(base, BASE)))
    stop("invalid datatype")
  if(pmatch(base, BASE) == -1)
    stop("ambiguous datatype")
  base <- match.arg(base, BASE)
  
  if(base=="size"){
    tmp <- invSize(x, datatype, size=level, conf=conf)
  }else if(base=="coverage"){
    tmp <- invChat(x, datatype, C=level, conf=conf)
  }
  
  tmp <- tmp[!duplicated(tmp),]
  
  nam <- names(x)
  if(is.null(nam)){
    tmp
  }else if(ncol(tmp)==6){
    tmp <- cbind(site=nam, tmp)
  }else{
    tmp <- cbind(site=rep(nam, each=3), tmp)
  }
  rownames(tmp) <- NULL
  tmp
}



# -----------------
# 2015-12-27, add transformation function 
# from incidence raw data to incidence frequencies data (iNEXT input format)
# 

###############################################
#' Transform incidence raw data to incidence frequencies (iNEXT input format) 
#' 
#' \code{as.incfreq}: transform incidence raw data (a species by sites presence-absence matrix) to incidence frequencies data (iNEXT input format, a row-sum frequencies vector contains total number of sampling units).
#' @param x a \code{data.frame} or \code{matirx} of species by sites presence-absence matrix.
#' @return a \code{vector} of species incidence frequencies, the first entry of the input data must be total number of sampling units.
#' @examples
#' data(ciliates)
#' lapply(ciliates, as.incfreq)
#' 
#' @export
#' 
as.incfreq <- function(x){
  class_x <- class(x)[1]
  if(class_x == "data.frame" | class_x == "matrix"){
    a <- sort(as.numeric(unique(c(unlist(x)))))
    if(!identical(a, c(0,1))){
      warning("Invalid data type, the element of species by sites presence-absence matrix should be 0 or 1. Set nonzero elements as 1.")
      x <- (x > 0)
    }
    nT <- ncol(x)
    y <- rowSums(x)
    y <- c(nT, y)
    # names(y) <- c("nT", rownames(x))
    y
  }else if(class_x=="numeric" | class_x=="integer" | class_x=="double"){
    warnings("Ambiguous data type, the input object is a vector. Set total number of sampling units as 1.")
    c(1, x) 
  }else{
    stop("Invalid data type, it should be a data.frame or matrix.")
  }
}

###############################################
#' Transform abundance raw data to abundance row-sum counts (iNEXT input format) 
#' 
#' \code{as.abucount}: transform abundance raw data (a species by sites matrix) to abundance rwo-sum counts data (iNEXT input format).
#' @param x a \code{data.frame} or \code{matirx} of species by sites matrix.
#' @return a \code{vector} of species abundance row-sum counts.
#' @examples
#' data(ciliates)
#' lapply(ciliates, as.abucount)
#' 
#' @export
#' 
as.abucount <- function(x){
  class_x <- class(x)[1]
  if(class_x == "data.frame" | class_x == "matrix"){
    y <- rowSums(x)
    y
  }else if(class_x=="numeric" | class_x=="integer" | class_x=="double"){
    warnings("Ambiguous data type, the input object is a vector. Set total number of sampling units as 1.")
    x 
  }else{
    stop("invalid data type, it should be a data.frame or matrix.")
  }
}