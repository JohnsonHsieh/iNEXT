estimateD <- function (x, q = c(0,1,2), datatype = "abundance", base = "size", level = NULL,nboot=50,
                          conf = NULL) 
{
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if (is.na(pmatch(datatype, TYPE))) 
    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1) 
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  if (datatype == "incidence") {
    stop("datatype=\"incidence\" was no longer supported after v2.0.8, \n         please try datatype=\"incidence_freq\".")
  }
  if (datatype == "incidence_freq") 
    datatype <- "incidence"
  if (datatype == "incidence_raw") {
    if (class(x) == "data.frame" | class(x) == "matrix") 
      x <- as.incfreq(x)
    else if (class(x) == "list") 
      x <- lapply(x, as.incfreq)
    datatype <- "incidence"
  }
  BASE <- c("size", "coverage")
  if (is.na(pmatch(base, BASE))) 
    stop("invalid datatype")
  if (pmatch(base, BASE) == -1) 
    stop("ambiguous datatype")
  base <- match.arg(base, BASE)
  if (base == "size") {
    tmp <- invSize_qs(x, q, datatype, size = level, nboot, conf = conf)
  }
  else if (base == "coverage") {
    tmp <- invChat_qs(x, q, datatype, C = level, nboot, conf = conf)
  }
  tmp <- tmp[!duplicated(tmp), ]
  tmp
}
invChat <- function (x, q, datatype = "abundance", C = NULL,nboot=50, conf = NULL) 
{
  TYPE <- c("abundance", "incidence")
  if (is.na(pmatch(datatype, TYPE))) 
    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1) 
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  if (class(x) == "numeric" | class(x) == "integer"){
    x <- list(data = x)
  }
  if (class(x) == "data.frame" | class(x) ==  "matrix"){
    datalist <- lapply(1:ncol(x), function(i) x[,i])
    if(is.null(colnames(x))) names(datalist) <-  paste0("data",1:ncol(x)) else names(datalist) <- colnames(x)
    x <- datalist
  }
  if (datatype == "abundance") {
    if (class(x) == "list") {
      if (is.null(C)) {
        C <- min(unlist(lapply(x, function(x) iChat.Ind(x,2*sum(x)))))
      }
      Community = rep(names(x),each = length(q))
      out <- do.call(rbind, lapply(x, function(x) invChat.Ind_qs(x, q, C,nboot, conf)))
      out$Community <- Community
      out <- out[,c(ncol(out),seq(1,(ncol(out)-1)))]
      rownames(out) <- NULL
    }else {
      stop("Wrong data format, dataframe/matrix or list would be accepted")
    }
  }else if (datatype == "incidence") {
    if (class(x) == "list") {
      if (is.null(C)) {
        C <- min(unlist(lapply(x, function(x) Chat.Sam(x,2*sum(x)))))
      }
      Community = rep(names(x),each = length(q))
      out <- do.call(rbind, lapply(x, function(x) invChat.Sam_qs(x,q,C,nboot, conf)))
      out$Community <- Community
      out <- out[,c(ncol(out),seq(1,(ncol(out)-1)))]
      rownames(out) <- NULL
    }else {
      stop("Wrong data format, dataframe/matrix or list would be accepted")
    }
  }
  out
}

invChat.Ind <- function (x, q, C, nboot=50, conf = NULL) 
{
  x <- x[x>0] ####added by yhc
  m <- NULL
  n <- sum(x)
  refC <- Chat.Ind(x, n)
  f <- function(m, C) abs(Chat.Ind(x, m) - C)
  if (refC > C) {
    opt <- optimize(f, C = C, lower = 0, upper = sum(x))
    mm <- opt$minimum
    mm <- round(mm)
  }
  if (refC <= C) {
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    if (f1 > 0 & f2 > 0) {
      A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
    }
    if (f1 > 1 & f2 == 0) {
      A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
    }
    if (f1 == 1 & f2 == 0) {
      A <- 1
    }
    if (f1 == 0 & f2 == 0) {
      A <- 1
    }
    mm <- (log(n/f1) + log(1 - C))/log(A) - 1
    mm <- n + mm
    mm <- round(mm)
  }
  if (mm > 2 * n) 
    warning("The maximum size of the extrapolation exceeds double reference sample size, the results for q = 0 may be subject to large prediction bias.")
  method <- ifelse(mm < n, "interpolated", ifelse(mm == 
                                                    n, "observed", "extrapolated"))
  if (nboot==0) {
    # SamCov = Chat.Ind(x, mm)
    # qD_coverage <- TD.m.est(x,mm,q)
    # out <- tibble(m=mm, Coverage = SamCov, q = q, Diversity = qD_coverage)
    out <- subset(iNEXT.Ind(x,q,m = c(1,mm),se = FALSE),m==mm)
    out <- out[,c(1,2,3,5,4)]
  }else {
    out <- subset(iNEXT.Ind(x,q,m = c(1,mm),se = TRUE,conf = conf,nboot = nboot), m==mm)
    out <- out[,c(1, 2, 3, 7, 4, 5, 6)]
  }
  out
}
invChat.Sam <- function (x, q, C, nboot=50, conf = NULL) 
{
  x <- x[x>0] ####added by yhc
  m <- NULL
  n <- max(x)
  refC <- Chat.Sam(x, n)
  f <- function(m, C) abs(Chat.Sam(x, m) - C)
  if (refC > C) {
    opt <- optimize(f, C = C, lower = 0, upper = max(x))
    mm <- opt$minimum
    mm <- round(mm)
  }
  if (refC <= C) {
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    U <- sum(x) - max(x)
    if (f1 > 0 & f2 > 0) {
      A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
    }
    if (f1 > 1 & f2 == 0) {
      A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
    }
    if (f1 == 1 & f2 == 0) {
      A <- 1
    }
    if (f1 == 0 & f2 == 0) {
      A <- 1
    }
    mm <- (log(U/f1) + log(1 - C))/log(A) - 1
    mm <- n + mm
    mm <- round(mm)
  }
  if (mm > 2 * n) 
    warning("The maximum size of the extrapolation exceeds double reference sample size, the results for q = 0 may be subject to large prediction bias.")
  if (is.null(conf)) {
    # SamCov = Chat.Sam(x, mm)
    # qD_coverage <- TD.m.est_inc(x,mm,q)
    # out <- tibble(t_=mm, Coverage = SamCov, q = q, Diversity = qD_coverage)
    out <- subset(iNEXT.Sam(x,q,t = c(1,mm),se = FALSE), t==mm)
    out <- out[,c(1,2,3,5,4)]
  }else {
    out <- subset(iNEXT.Sam(x,q,t = c(1,mm),se = TRUE,conf = conf,nboot = nboot), t==mm)
    out <- out[, c(1, 2, 3, 7, 4, 5, 6)]
  }
  out
}

invSize <- function(x, q, datatype="abundance", size=NULL, nboot=50, conf=NULL){
  TYPE <- c("abundance", "incidence")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  class_x <- class(x)[1]
  if (class(x) == "numeric" | class(x) == "integer"){
    x <- list(data = x)
  }
  if (class(x) == "data.frame" | class(x) ==  "matrix"){
    datalist <- lapply(1:ncol(x), function(i) x[,i])
    if(is.null(colnames(x))) names(datalist) <-  paste0("data",1:ncol(x)) else names(datalist) <- colnames(x)
    x <- datalist
  }
  if(datatype=="abundance"){
    if (class(x) == "list") {
      if (is.null(size)) {
        size <- min(unlist(lapply(x, function(x) sum(x))))
      } 
      Community = rep(names(x),each = length(q))
      out <- do.call(rbind, lapply(x, function(x) invSize.Ind_qs(x,q,size,nboot,conf)))
      out$Community <- Community
      out <- out[,c(ncol(out),seq(1,(ncol(out)-1)))]
      rownames(out) <- NULL
    }else {
      stop("Wrong data format, dataframe/matrix or list would be accepted")
    }
  }else if (datatype == "incidence") {
    if (class(x) == "list") {
      if (is.null(size)) {
        size <- min(unlist(lapply(x, function(x) max(x))))
      }
      Community = rep(names(x),each = length(q))
      out <- do.call(rbind, lapply(x, function(x) invSize.Sam_qs(x,q,size,nboot,conf)))
      out$Community <- Community
      out <- out[,c(ncol(out),seq(1,(ncol(out)-1)))]
      rownames(out) <- NULL
    }else {
      stop("Wrong data format, dataframe/matrix or list would be accepted")
    }
  }
  out
}


invSize.Ind <- function(x, q, size, nboot=50, conf=NULL){
  m <- NULL # no visible binding for global variable 'm'
  
  if(is.null(size)){
    size <- sum(x)
  }
  if(nboot==0){
    method <- ifelse(size<sum(x), "interpolated", ifelse(size==sum(x), "observed", "extrapolated"))
    out <- subset(iNEXT.Ind(x,q,m = c(1,size),se = FALSE), m==size)
    out <- out[,c(1,2,3,5,4)]
    # out <- data.frame(m=size, method=method, 
    #                   SamCov=round(Chat.Ind(x,size),3),
    #                   SpeRic=round(Dqhat.Ind(x,0,size),3),
    #                   ShaDiv=round(Dqhat.Ind(x,1,size),3),
    #                   SimDiv=round(Dqhat.Ind(x,2,size),3))
    # colnames(out) <- c("m", "method", "SC", "q = 0", "q = 1", "q = 2")
  }else{
    out <- subset(iNEXT.Ind(x,q,m = c(1,size),se = TRUE,conf = conf,nboot = nboot), m==size)
    out <- out[,c(1, 2, 3, 7, 4, 5, 6)]
  }
  out
}

invSize.Sam <- function(x, q, size, nboot=50, conf=NULL){
  m <- NULL # no visible binding for global variable 'm'
  
  if(is.null(size)){
    size <- max(x)
  }
  if(nboot==0){
    method <- ifelse(size<max(x), "interpolated", ifelse(size==max(x), "observed", "extrapolated"))
    out <- subset(iNEXT.Sam(x,q,t = c(1,size),se = FALSE), t==size)
    out <- out[,c(1,2,3,5,4)]
    # out <- data.frame(t=size, method=method, 
    #                   SamCov=round(Chat.Sam(x,size),3),
    #                   SpeRic=round(Dqhat.Sam(x,0,size),3),
    #                   ShaDiv=round(Dqhat.Sam(x,1,size),3),
    #                   SimDiv=round(Dqhat.Sam(x,2,size),3))
    # colnames(out) <- c("m", "method", "SC", "q = 0", "q = 1", "q = 2")
  }else{
    out <- subset(iNEXT.Sam(x,q,t = c(1,size),se = TRUE,conf = conf,nboot = nboot), t==size)
    out <- out[, c(1, 2, 3, 7, 4, 5, 6)]
    # tmp0 <- iNEXT.Sam(x, q=0, t=c(1,size), se = TRUE, conf = conf)
    # tmp1 <- iNEXT.Sam(x, q=1, t=c(1,size), se = TRUE, conf = conf)
    # tmp2 <- iNEXT.Sam(x, q=2, t=c(1,size), se = TRUE, conf = conf)
    # tmp <- subset(rbind(tmp0,tmp1,tmp2), t==size)
    # out <- tmp[,c(1,2,3,7,4,5,6)]
    # out[,4:7] <- round(out[,4:7],3)
  }
  out
}
