invChat.Ind <- function(x, C, upper=sum(x)){
  f <- function(m, C) abs(Chat.Ind(x,m)-C)
  opt <- optimize(f, C=C, lower=0, upper=2*sum(x))
  mm <- opt$minimum
  if(opt$objective>0.0001){
    warning("The specific extrapolation is too far away to reference sample size, the return won't be robust.")
    m <- 1:(2*sum(x))
    mm <- predict(smooth.spline(x=Chat.Ind(x,m), y=m), x=C)$y
  }
  n <- sum(x)
  method <- ifelse(mm<n, "rarefied", ifelse(mm==n, "observed", "extrapolation"))
  mm <- ceiling(mm)
  out <- data.frame(m=mm, method=method, 
                    SamCov=round(Chat.Ind(x,mm),3), 
                    SpeRic=round(Dqhat.Ind(x,0,mm),3),
                    ShaDiv=round(Dqhat.Ind(x,1,mm),3),
                    SimDiv=round(Dqhat.Ind(x,2,mm),3))
  out  
}


invChat.Sam <- function(x, C, upper=max(x)){
  f <- function(m, C) abs(Chat.Sam(x,m)-C)
  opt <- optimize(f, C=C, lower=0, upper=2*sum(x))
  mm <- opt$minimum
  if(opt$objective>0.0001){
    warning("The specific extrapolation is too far away to reference sample size, the return won't be robust.")
    m <- 1:(2*sum(x))
    mm <- predict(smooth.spline(x=Chat.Sam(x,m), y=m), x=C)$y
  }
  n <- max(x)
  method <- ifelse(mm<n, "rarefied", ifelse(mm==n, "observed", "extrapolation"))
  mm <- ceiling(mm)
  out <- data.frame(t=mm, method=method, 
                    SamCov=round(Chat.Sam(x,mm),3),
                    SpeRic=round(Dqhat.Sam(x,0,mm),3),
                    ShaDiv=round(Dqhat.Sam(x,1,mm),3),
                    SimDiv=round(Dqhat.Sam(x,2,mm),3))
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
      method <- ifelse(size<sum(x), "rarefied", ifelse(size==sum(x), "observed", "extrapolation"))
      data.frame(m=size, method=method, 
                 SamCov=round(Chat.Ind(x,size),3),
                 SpeRic=round(Dqhat.Ind(x,0,size),3),
                 ShaDiv=round(Dqhat.Ind(x,1,size),3),
                 SimDiv=round(Dqhat.Ind(x,2,size),3))
    }
    else if(class(x)=="list"){
      if(is.null(size)){
        size <- min(unlist(lapply(x, function(x) sum(x))))
      }    
      
      do.call(rbind, lapply(x, function(x){
        method <- ifelse(size<sum(x), "rarefied", ifelse(size==sum(x), "observed", "extrapolation"))
        data.frame(m=size, method=method, 
                   SamCov=round(Chat.Ind(x,size),3),
                   SpeRic=round(Dqhat.Ind(x,0,size),3),
                   ShaDiv=round(Dqhat.Ind(x,1,size),3),
                   SimDiv=round(Dqhat.Ind(x,2,size),3))
      }))
    }else if(class(x)=="data.fram" | class(x)=="matrix"){
      if(is.null(size)){
        size <- min(unlist(apply(x, 2, function(x) sum(x))))
      }
      do.call(rbind, apply(x, 2, function(x){
        method <- ifelse(size<sum(x), "rarefied", ifelse(size==sum(x), "observed", "extrapolation"))
        data.frame(m=size, method=method, 
                   SamCov=round(Chat.Ind(x,size),3),
                   SpeRic=round(Dqhat.Ind(x,0,size),3),
                   ShaDiv=round(Dqhat.Ind(x,1,size),3),
                   SimDiv=round(Dqhat.Ind(x,2,size),3))
      }))
    }
  }else if(datatype=="incidence"){
    
    if(class(x)=="numeric" | class(x)=="integer"){
      if(is.null(size)){
        size <- max(x)
      }
      method <- ifelse(size<max(x), "rarefied", ifelse(size==max(x), "observed", "extrapolation"))
      data.frame(t=size, method=method, 
                 SamCov=round(Chat.Sam(x,size),3),
                 SpeRic=round(Dqhat.Sam(x,0,size),3),
                 ShaDiv=round(Dqhat.Sam(x,1,size),3),
                 SimDiv=round(Dqhat.Sam(x,2,size),3))
    }
    else if(class(x)=="list"){
      if(is.null(size)){
        size <- min(unlist(lapply(x, function(x) max(x))))
      }
      do.call(rbind, lapply(x, function(x){
        method <- ifelse(size<max(x), "rarefied", ifelse(size==max(x), "observed", "extrapolation"))
        data.frame(t=size, method=method, 
                   SamCov=round(Chat.Sam(x,size),3),
                   SpeRic=round(Dqhat.Sam(x,0,size),3),
                   ShaDiv=round(Dqhat.Sam(x,1,size),3),
                   SimDiv=round(Dqhat.Sam(x,2,size),3))
      }))
    }else if(class(x)=="data.fram" | class(x)=="matrix"){
      if(is.null(size)){
        size <- min(unlist(apply(x, 2, function(x) max(x))))
      }
      do.call(rbind, apply(x, 2, function(x){
        method <- ifelse(size<max(x), "rarefied", ifelse(size==max(x), "observed", "extrapolation"))
        data.frame(t=size, method=method, 
                   SamCov=round(Chat.Sam(x,size),3),
                   SpeRic=round(Dqhat.Sam(x,0,size),3),
                   ShaDiv=round(Dqhat.Sam(x,1,size),3),
                   SimDiv=round(Dqhat.Sam(x,2,size),3))
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
#' Compute species diversity with fixed sample coverage
#' 
#' \code{Auxiliary} compute species diversity with fixed sample size or sample coverage.
#' @param x a \code{data.frame} or \code{list} for species abundance or incidence frequencies.
#' @param datatype the data type of input data. That is individual-based abundance data (\code{datatype = "abundance"}) or sample-based incidence data (\code{datatype = "incidence"}).
#' @param base the compared base, That is sample-size-based (\code{base="size"}) or coverage-based (\code{base="coverage"}).
#' @param ref the compared reference sample size, sampling units or sample coverage. 
#' If \code{base="size"} and \code{red=NULL}, this function compute the minimum level of sample size or samplping units for all sites. 
#' If \code{base="coverage"} and \code{ref=NULL}, this function compute the minimum level of sample coverage. 
#' @return a \code{data.frame} of species diversity table with fixed minimum reference sample size or minimum sample coverage for all sites.
#' @examples
#' data(spider)
#' Auxiliary(spider, "abundance", base="size")
#' Auxiliary(spider, "abundance", base="coverage")
#' 
#' data(ant)
#' Auxiliary(ant, "incidence", base="coverage", ref=0.985)
#' @export
Auxiliary <- function(x, datatype="abundance", base="size", ref=NULL){
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
    invSize(x, datatype, size=ref)
  }else if(base=="coverage"){
    invChat(x, datatype, C=ref)
  }
}

# Compare(spider, datatype="abundance", base="size")
# Compare(spider, datatype="abundance", base="coverage")
# Compare(ant, datatype="incidence", base="size")
# Compare(ant, datatype="incidence", base="coverage")
