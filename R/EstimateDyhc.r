source("iNEXT_generalq.r")
estimateDyhc <- function (x, q, datatype = "abundance", base = "size", level = NULL, 
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
    tmp <- invSize(x, datatype, size = level, conf = conf)
  }
  else if (base == "coverage") {
    tmp <- invChatyhc(x,q, datatype, C = level, conf = conf)
  }
  tmp <- tmp[!duplicated(tmp), ]
  tmp
}
invChatyhc <- function (x, q, datatype = "abundance", C = NULL, conf = NULL) 
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
        C <- min(unlist(lapply(x, function(x) iNEXT:::Chat.Ind(x,2*sum(x)))))
      }
      Community = rep(names(x),each = length(q))
      out <- do.call(rbind, lapply(x, function(x) invChat.Indyhc(x, q, C, conf))) %>%
        mutate(Community = Community)
    }else {
      stop("Wrong data format, dataframe/matrix or list would be accepted")
    }
  }else if (datatype == "incidence") {
    if (class(x) == "list") {
      if (is.null(C)) {
        C <- min(unlist(lapply(x, function(x) iNEXT:::Chat.Sam(x,2*sum(x)))))
      }
      Community = rep(names(x),each = length(q))
      out <- do.call(rbind, lapply(x, function(x) invChat.Samyhc(x,q,C, conf))) %>%
        mutate(Community = Community)
    }
  }
  out
}

invChat.Indyhc <- function (x, q, C, conf = NULL) 
{
  x <- x[x>0] ####added by yhc
  m <- NULL
  n <- sum(x)
  refC <- iNEXT:::Chat.Ind(x, n)
  f <- function(m, C) abs(iNEXT:::Chat.Ind(x, m) - C)
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
  if (is.null(conf)) {
    SamCov = iNEXT:::Chat.Ind(x, mm)
    qD_coverage <- TD.m.est_yhc(x,mm,q)
    out <- tibble(m=mm, Coverage = SamCov, q = q, Diversity = qD_coverage)
  }
  ###==== bootstrap is still under construction====
  # else {
  #   tmp0 <- iNEXT.Ind(x, q = 0, m = c(1, mm), se = TRUE, 
  #                     conf = conf)
  #   tmp1 <- iNEXT.Ind(x, q = 1, m = c(1, mm), se = TRUE, 
  #                     conf = conf)
  #   tmp2 <- iNEXT.Ind(x, q = 2, m = c(1, mm), se = TRUE, 
  #                     conf = conf)
  #   tmp <- subset(rbind(tmp0, tmp1, tmp2), m == mm)
  #   out <- tmp[, c(1, 2, 3, 7, 4, 5, 6)]
  #   out[, 4:7] <- round(out[, 4:7], 3)
  # }
  ###==== bootstrap is still under construction====
  out
}
invChat.Samyhc <- function (x, q, C, conf = NULL) 
{
  x <- x[x>0] ####added by yhc
  m <- NULL
  n <- max(x)
  refC <- iNEXT:::Chat.Sam(x, n)
  f <- function(m, C) abs(iNEXT:::Chat.Sam(x, m) - C)
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
    SamCov = iNEXT:::Chat.Sam(x, mm)
    qD_coverage <- TD.m.est_inc_yhc(x,mm,q)
    out <- tibble(t_=mm, Coverage = SamCov, q = q, Diversity = qD_coverage)
  }
  ###==== bootstrap is still under construction====
  # else {
  #   tmp0 <- iNEXT.Sam(x, q = 0, t = c(1, mm), se = TRUE, 
  #                     conf = conf)
  #   tmp1 <- iNEXT.Sam(x, q = 1, t = c(1, mm), se = TRUE, 
  #                     conf = conf)
  #   tmp2 <- iNEXT.Sam(x, q = 2, t = c(1, mm), se = TRUE, 
  #                     conf = conf)
  #   tmp <- subset(rbind(tmp0, tmp1, tmp2), t == mm)
  #   out <- tmp[, c(1, 2, 3, 7, 4, 5, 6)]
  #   out[, 4:7] <- round(out[, 4:7], 3)
  # }
  ###==== bootstrap is still under construction====
  out
}
