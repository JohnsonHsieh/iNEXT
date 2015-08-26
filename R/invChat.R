invChat.Ind <- function(x, C)
{
  n <- sum(x)
  refC <- Chat.Ind(x,n)
  f <- function(m, C) abs(Chat.Ind(x,m)-C)
  if(refC >= C)
  {
     opt <- optimize(f, C=C, lower=0, upper=sum(x))
     mm <- opt$minimum
     mm <- round(mm)
  }
  if(refC < C)
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
  method <- ifelse(mm<n, "racefaction", ifelse(mm==n, "observation", "extrapolation"))
  out <- data.frame(m=mm, method=method, 
                    "SC"=round(Chat.Ind(x,mm),4), 
                    "q = 0"=round(Dqhat.Ind(x,0,mm),3),
                    "q = 1"=round(Dqhat.Ind(x,1,mm),3),
                    "q = 2"=round(Dqhat.Ind(x,2,mm),3))
  colnames(out) <- c("m", "method", "SC", "q = 0", "q = 1", "q = 2")
  out
}

invChat.Sam <- function(x, C)
{
  n <- max(x)
  refC <- Chat.Sam(x,n)
  f <- function(m, C) abs(Chat.Sam(x,m)-C)
  if(refC >= C)
  {
     opt <- optimize(f, C=C, lower=0, upper=max(x))
     mm <- opt$minimum
     mm <- round(mm)
  }
  if(refC < C)
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
  
  method <- ifelse(mm<n, "rarefaction", ifelse(mm==n, "observation", "extrapolation"))
  out <- data.frame(t=mm, method=method, 
                    "SC"=round(Chat.Ind(x,mm),4), 
                    "q = 0"=round(Dqhat.Ind(x,0,mm),3),
                    "q = 1"=round(Dqhat.Ind(x,1,mm),3),
                    "q = 2"=round(Dqhat.Ind(x,2,mm),3))
  colnames(out) <- c("t", "method", "SC", "q = 0", "q = 1", "q = 2")
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

invChat <- function(x, datatype="abundance", C=NULL){
  TYPE <- c("abundance", "incidence")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  
  if(datatype=="abundance"){
    if(class(x)=="numeric" | class(x)=="integer"){
      if(is.null(C)){
        C <- Chat.Ind(x, 2*sum(x))
      }
      invChat.Ind(x, C)
    }
    else if(class(x)=="list"){
      if(is.null(C)){
        C <- min(unlist(lapply(x, function(x) Chat.Ind(x, sum(x)))))
      }
      do.call(rbind, lapply(x, function(x) invChat.Ind(x, C)))
    }else if(class(x)=="data.fram" | class(x)=="matrix"){
      if(is.null(C)){
        C <- min(unlist(apply(x, 2, function(x) Chat.Ind(x, sum(x)))))
      }
      do.call(rbind, apply(x, 2, function(x) invChat.Ind(x, C)))
    }
  }else if(datatype=="incidence"){
    if(class(x)=="numeric" | class(x)=="integer"){
      if(is.null(C)){
        C <- Chat.Sam(x, 2*max(x))
      }
      invChat.Ind(x, C)
    }
    else if(class(x)=="list"){
      if(is.null(C)){
        C <- min(unlist(lapply(x, function(x) Chat.Sam(x, max(x)))))
      }
      do.call(rbind, lapply(x, function(x) invChat.Sam(x, C)))
    }else if(class(x)=="data.fram" | class(x)=="matrix"){
      if(is.null(C)){
        C <- min(unlist(apply(x, 2, function(x) Chat.Sam(x, max(x)))))
      }
      do.call(rbind, apply(x, 2, function(x) invChat.Sam(x, C)))
    }
    
  }
}

invSize <- function(x, datatype="abundance", size=NULL){
  TYPE <- c("abundance", "incidence")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  if(datatype=="abundance"){
    if(class(x)=="numeric" | class(x)=="integer"){
      if(is.null(size)){
        size <- sum(x)
      }
      method <- ifelse(size<sum(x), "racefaction", ifelse(size==sum(x), "observation", "extrapolation"))
      out <- data.frame(m=size, method=method, 
                 SamCov=round(Chat.Ind(x,size),3),
                 SpeRic=round(Dqhat.Ind(x,0,size),3),
                 ShaDiv=round(Dqhat.Ind(x,1,size),3),
                 SimDiv=round(Dqhat.Ind(x,2,size),3))
	  colnames(out) <- c("m", "method", "SC", "q = 0", "q = 1", "q = 2")
	  out
    }
    else if(class(x)=="list"){
      if(is.null(size)){
        size <- min(unlist(lapply(x, function(x) sum(x))))
      }    
      
      do.call(rbind, lapply(x, function(x){
        method <- ifelse(size<sum(x), "racefaction", ifelse(size==sum(x), "observation", "extrapolation"))
        out <- data.frame(m=size, method=method, 
                   SamCov=round(Chat.Ind(x,size),3),
                   SpeRic=round(Dqhat.Ind(x,0,size),3),
                   ShaDiv=round(Dqhat.Ind(x,1,size),3),
                   SimDiv=round(Dqhat.Ind(x,2,size),3))
		colnames(out) <- c("m", "method", "SC", "q=0", "q=1", "q=2")
		out
      }))
    }else if(class(x)=="data.fram" | class(x)=="matrix"){
      if(is.null(size)){
        size <- min(unlist(apply(x, 2, function(x) sum(x))))
      }
      do.call(rbind, apply(x, 2, function(x){
        method <- ifelse(size<sum(x), "racefaction", ifelse(size==sum(x), "observation", "extrapolation"))
        out <- data.frame(m=size, method=method, 
                   SamCov=round(Chat.Ind(x,size),3),
                   SpeRic=round(Dqhat.Ind(x,0,size),3),
                   ShaDiv=round(Dqhat.Ind(x,1,size),3),
                   SimDiv=round(Dqhat.Ind(x,2,size),3))
		colnames(out) <- c("m", "method", "SC", "q=0", "q=1", "q=2")
		out
      }))
    }
  }else if(datatype=="incidence"){
    
    if(class(x)=="numeric" | class(x)=="integer"){
      if(is.null(size)){
        size <- max(x)
      }
      method <- ifelse(size<max(x), "racefaction", ifelse(size==max(x), "observation", "extrapolation"))
      out <- data.frame(t=size, method=method, 
                 SamCov=round(Chat.Sam(x,size),3),
                 SpeRic=round(Dqhat.Sam(x,0,size),3),
                 ShaDiv=round(Dqhat.Sam(x,1,size),3),
                 SimDiv=round(Dqhat.Sam(x,2,size),3))
	  colnames(out) <- c("m", "method", "SC", "q=0", "q=1", "q=2")
	  out
    }
    else if(class(x)=="list"){
      if(is.null(size)){
        size <- min(unlist(lapply(x, function(x) max(x))))
      }
      do.call(rbind, lapply(x, function(x){
        method <- ifelse(size<max(x), "racefaction", ifelse(size==max(x), "observation", "extrapolation"))
        out <- data.frame(t=size, method=method, 
                   SamCov=round(Chat.Sam(x,size),3),
                   SpeRic=round(Dqhat.Sam(x,0,size),3),
                   ShaDiv=round(Dqhat.Sam(x,1,size),3),
                   SimDiv=round(Dqhat.Sam(x,2,size),3))
		colnames(out) <- c("m", "method", "SC", "q=0", "q=1", "q=2")
		out
      }))
    }else if(class(x)=="data.fram" | class(x)=="matrix"){
      if(is.null(size)){
        size <- min(unlist(apply(x, 2, function(x) max(x))))
      }
      do.call(rbind, apply(x, 2, function(x){
        method <- ifelse(size<max(x), "racefaction", ifelse(size==max(x), "observation", "extrapolation"))
        out <- data.frame(t=size, method=method, 
                   SamCov=round(Chat.Sam(x,size),3),
                   SpeRic=round(Dqhat.Sam(x,0,size),3),
                   ShaDiv=round(Dqhat.Sam(x,1,size),3),
                   SimDiv=round(Dqhat.Sam(x,2,size),3))
		colnames(out) <- c("m", "method", "SC", "q=0", "q=1", "q=2")
		out
      }))
    }
    
  }
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
#' Compute species diversity with specific level
#' 
#' \code{Auxiliary}: computes species diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage.
#' @param x a \code{data.frame} or \code{list} of species abundances or incidence frequencies.\cr 
#' If \code{datatype = "incidence"}, then the first entry of the input data must be total number of sampling units, followed 
#' by species incidence frequencies in each column or list.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or 
#' sampling-unit-based incidence data (\code{datatype = "incidence"}).
#' @param base comparison base: sample-size-based (\code{base="size"}) or coverage-based \cr (\code{base="coverage"}).
#' @param level an integer specifying a particular sample size or a number (between 0 and 1) specifying a particular value of sample coverage. 
#' If \code{base="size"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample size among all sites. 
#' If \code{base="coverage"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample coverage among all sites. 
#' @return a \code{data.frame} of species diversity table including the reference sample size, sample coverage,
#' method (rarefaction or extrapolation), and diversity estimates with q = 0, 1, and 2 for the user-specified sample size or sample coverage.
#' @examples
#' data(spider)
#' Auxiliary(spider, "abundance", base="size", level=NULL)
#' Auxiliary(spider, "abundance", base="coverage", level=NULL)
#' 
#' data(ant)
#' Auxiliary(ant, "incidence", base="coverage", level=0.985)
#' @export
Auxiliary <- function(x, datatype="abundance", base="size", level=NULL){
  TYPE <- c("abundance", "incidence")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  
  BASE <- c("size", "coverage")
  if(is.na(pmatch(base, BASE)))
    stop("invalid datatype")
  if(pmatch(base, BASE) == -1)
    stop("ambiguous datatype")
  base <- match.arg(base, BASE)
  
  if(base=="size"){
    invSize(x, datatype, size=level)
  }else if(base=="coverage"){
    invChat(x, datatype, C=level)
  }
}

# Compare(spider, datatype="abundance", base="size")
# Compare(spider, datatype="abundance", base="coverage")
# Compare(ant, datatype="incidence", base="size")
# Compare(ant, datatype="incidence", base="coverage")
