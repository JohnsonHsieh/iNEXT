#
#
###########################################
#' Exhibit basic data information
#' 
#' \code{DataInfo}: exhibits basic data information
#' 
#' @param x a vector/matrix/data.frame/list of species abundances or incidence frequencies.\cr If \code{datatype = "incidence_freq"}, 
#' then the first entry of the input data must be total number of sampling units, followed by species incidence frequencies.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @return a data.frame of basic data information including sample size, observed species richness, sample coverage estimate, and the first ten abundance/incidence frequency counts.
#' @examples 
#' data(spider)
#' DataInfo(spider, datatype="abundance")
#' @export
DataInfo <- function(x, datatype="abundance"){
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  
  if(datatype=="incidence_freq") datatype <- "incidence"
  
  if(datatype=="incidence_raw"){
    if(inherits(x, "list")){
      x <- lapply(x, as.incfreq)
    }else{
      x <- as.incfreq(x)
    }
    datatype <- "incidence"
  }
  
  Fun.abun <- function(x){
    n <- sum(x)
    fk <- sapply(1:10, function(k) sum(x==k))
    f1 <- fk[1]
    f2 <- fk[2]
    Sobs <- sum(x>0)
    f0.hat <- ifelse(f2==0, (n-1)/n*f1*(f1-1)/2, (n-1)/n*f1^2/2/f2)  #estimation of unseen species via Chao1
    A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
    Chat <- round(1 - f1/n*A, 4)
    c(n, Sobs, Chat, fk)
  }
  
  Fun.ince <- function(x){
    nT <- x[1]
    x <- x[-1]
    U <- sum(x)
    Qk <- sapply(1:10, function(k) sum(x==k))
    Q1 <- Qk[1]
    Q2 <- Qk[2]
    Sobs <- sum(x>0)
    Q0.hat <- ifelse(Q2==0, (nT-1)/nT*Q1*(Q1-1)/2, (nT-1)/nT*Q1^2/2/Q2)  #estimation of unseen species via Chao2
    A <- ifelse(Q1>0, nT*Q0.hat/(nT*Q0.hat+Q1), 1)
    Chat <- round(1 - Q1/U*A,4)
    out <- c(nT, U, Sobs, Chat, Qk)
  }
  
  if(datatype == "abundance"){
    if(inherits(x, "numeric") | inherits(x, "integer")){
      out <- matrix(Fun.abun(x), nrow=1)
    }else if(inherits(x, "list")){
      out <- do.call("rbind", lapply(x, Fun.abun))
    } else if(inherits(x, "matrix") | inherits(x, "data.frame")){
      out <- t(apply(as.matrix(x), 2, Fun.abun))  
    }
    if(nrow(out) > 1 | inherits(x, "list")){
      out <- data.frame(site=rownames(out), out)
      colnames(out) <-  c("Assemblage", "n", "S.obs", "SC", paste("f",1:10, sep=""))
      rownames(out) <- NULL
    }else{
      out <- data.frame(site="site.1", out)
      colnames(out) <-  c("Assemblage", "n", "S.obs", "SC", paste("f",1:10, sep=""))
    }
    as.data.frame(out)
  }else if(datatype == "incidence"){
    if(inherits(x, "numeric") | inherits(x, "integer")){
      out <- matrix(Fun.ince(x), nrow=1)
    }else if(inherits(x, "list")){
      out <- do.call("rbind", lapply(x, Fun.ince))
    } else if(inherits(x, "matrix") | inherits(x, "data.frame")){
      out <- t(apply(as.matrix(x), 2, Fun.ince))  
    }
    if(nrow(out) > 1 | inherits(x, "list")){
      out <- data.frame(site=rownames(out), out)
      colnames(out) <-  c("Assemblage","T", "U", "S.obs", "SC", paste("Q",1:10, sep=""))
      rownames(out) <- NULL
    }else{
      out <- data.frame(site="site.1", out)
      colnames(out) <-  c("Assemblage","T", "U", "S.obs", "SC", paste("Q",1:10, sep=""))
    }
    as.data.frame(out)
  }
}




#
#
###########################################
#' Estimation of species richness
#' 
#' \code{ChaoRichness}: estimation of species richness based on the methods proposed in Chao (1984, 1987)
#' 
#' @param x a \code{matrix}, \code{data.frame} (species by sites), or \code{list} of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, 
#' then the first entry of the input data must be total number of sampling units, followed by species incidence frequencies. 
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param conf a positive number \eqn{\le} 1 specifying the level of confidence interval. 
#' @return A data.frame of observed species richness, species richness estimate, s.e. and the associated confidence interval.
#' @seealso \code{\link{ChaoShannon}, \link{ChaoSimpson}}
#' @examples 
#' data(spider)
#' ChaoRichness(spider$Girdled, datatype="abundance")
#' @references 
#' Chao, A. (1984) Nonparametric estimation of the number of classes in a population. Scandinavian Journal of Statistics, 11, 265-270.\cr\cr
#' Chao, A. (1987) Estimating the population size for capture-recapture data with unequal catchability. Biometrics, 43, 783-791.
#' 
#' @export

ChaoRichness=function(x, datatype="abundance", conf=0.95){  
  
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid data type")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous data type")
  datatype <- match.arg(datatype, TYPE)
  
  if(datatype=="incidence_freq") datatype <- "incidence"
  
  if(datatype=="incidence_raw"){
    if(inherits(x, "list")){
      x <- lapply(x, as.incfreq)
    }else{
      x <- as.incfreq(x)
    }
    datatype <- "incidence"
  }
  
  
  if (!is.numeric(conf) || conf > 1 || conf < 0) {
    warning("\"conf\"(confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!")
    conf <- 0.95
  }
  
  myFun <- function(x){
    if(!is.numeric(x) & !is.matrix(x) & !is.data.frame(x))
      stop("invalid data structure")
    if(is.matrix(x) | is.data.frame(x)){
      if (ncol(x) != 1 & nrow(x) != 1)
        stop("invalid data structure")
    }
    z <- qnorm(1 - (1 - conf)/2)
    if(datatype=="abundance"){
      n <- sum(x)
      D <- sum(x>0)
    }else if(datatype=="incidence"){
      n <- x[1]
      x <- x[-1]
      D <- sum(x>0)
      U <- sum(x)    
    }
    
    f1=sum(x==1)
    f2=sum(x==2)
    if (f1 > 0 & f2 > 0)
    {
      S_Chao1 <- D + (n - 1)/n*f1^2/(2*f2)
      var_Chao1 <- f2*((n - 1)/n*(f1/f2)^2/2 + 
                         ((n - 1)/n)^2*(f1/f2)^3 + ((n - 1)/n)^2*(f1/f2)^4/4)
      
      t <- S_Chao1 - D
      K <- exp(z*sqrt(log(1 + var_Chao1/t^2)))
      CI_Chao1 <- c(D + t/K, D + t*K)
    } 
    else if (f1 > 1 & f2 == 0)
    {
      S_Chao1 <- D + (n - 1)/n*f1*(f1 - 1)/(2*(f2 + 1))
      var_Chao1 <- (n - 1)/n*f1*(f1 - 1)/2 + 
        ((n - 1)/n)^2*f1*(2*f1 - 1)^2/4 - ((n - 1)/n)^2*f1^4/4/S_Chao1
      
      t <- S_Chao1 - D
      K <- exp(z*sqrt(log(1 + var_Chao1/t^2)))
      CI_Chao1 <- c(D + t/K, D + t*K)
    } 
    else 
    {
      S_Chao1 <- D
      i <- c(1:max(x))
      i <- i[unique(x)]
      var_obs <- sum(sapply(i, function(i) sum(x==i)*(exp(-i) - exp(-2*i)))) - 
        (sum(sapply(i, function(i)i*exp(-i)*sum(x==i))))^2/n
      var_Chao1 <- var_obs
      P <- sum(sapply(i, function(i) sum(x==i)*exp(-i)/D))
      CI_Chao1 <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
    }
    out <- c(D, S_Chao1, var_Chao1^0.5, CI_Chao1[1], CI_Chao1[2])
    out[-1] <- round(out[-1],3)
    out <- data.frame(t(out))
    colnames(out) <- c("Observed", "Estimator", "Est_s.e.", paste(conf*100,"% Lower",sep=""), paste(conf*100,"% Upper",sep=""))
    return(out)
  }
  
  if(inherits(x, "numeric")){
    out <- myFun(x)
  }else if(inherits(x, "integer")){
	out <- myFun(x)
  }else if(inherits(x, "list")){
    out <- do.call("rbind", lapply(x, myFun))
  } else if(inherits(x, "matrix") | inherits(x, "data.frame")){
    out <- do.call("rbind", apply(as.matrix(x), 2, myFun))  
  }
  return(out)
}


#
#
###########################################
# Estimation of species richness
# 
# \code{BootstrapFun} Estimation of species richness via Chao (1984, 1987)
# 
# @param x a vector of species frequencies.
# @param FunName the R function to estimate the traget index.
# @param datatype the data type of input data. That is individual-based abundance data (\code{datatype = "abundance"}) or sample-based incidence data (\code{datatype = "incidence"}).
# @param B the number of bootstrap resampling times, default is \code{200}.
# @return standard error of the estimator
# @author Y.H. Lee
BootstrapFun <- function(x, FunName, datatype, B){
  
  if(!is.numeric(x) & !is.matrix(x) & !is.data.frame(x))
    stop("invalid data structure")
  if(is.matrix(x) | is.data.frame(x)){
    if (ncol(x) != 1 & nrow(x) != 1)
      stop("invalid data structure")
  }
  
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid data type")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous data type")
  datatype <- match.arg(datatype, TYPE)
  
  if(datatype=="incidence_freq") datatype <- "incidence"
  
  if(datatype=="incidence_raw"){
    if(inherits(x, "list")){
      x <- lapply(x, as.incfreq)
    }else{
      x <- as.incfreq(x)
    }
    datatype <- "incidence"
  }
  
  BootstrapFun.abun <- function(x, FunName, datatype, B) {
    n <- sum(x)
    f1 <- sum(x==1)
    f2 <- sum(x==2)
    f0.hat <- ifelse(f2==0, (n-1)/n*f1*(f1-1)/2, (n-1)/n*f1^2/2/f2)  #estimation of unseen species via Chao1
    A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
    Chat <- 1 - f1/n*A
    f0 <- max(round(f0.hat), 1)
    
    if(f0.hat==0){
      lambda <- 0
      if(sum(x>0)==1){
        warning("The Bootstrap community has only one species. Estimation is not robust.")
      }
    }else{
      lambda <- (1 - Chat) / sum(x / n * (1 - x / n)^n)
    }
    
    pi <- x / n * (1 - lambda * (1 - x /n)^n)
    pi.star <- c(pi, rep((1 - Chat) / f0, f0))
    #   set.seed(1234)
    X <- rmultinom(B, n, pi.star)
    se <- sd(apply(X, 2, function(x) FunName(x, datatype)))
    return(se)
  }
  
  BootstrapFun.ince <- function(y, FunName, datatype, B) {
    t <- y[1]
    y <- y[-1]
    y <- y[y>0]
    U <- sum(y)
    Q1 <- sum(y==1)
    Q2 <- sum(y==2)
    Q0.hat <- ifelse(Q2==0, (t-1)/t*Q1*(Q1-1)/2, (t-1)/t*Q1^2/2/Q2)  #estimation of unseen species via Chao2
    A <- ifelse(Q1>0, t*Q0.hat/(t*Q0.hat+Q1), 1)
    Chat <- 1 - Q1/U*A
    Q0 <- max(round(Q0.hat), 1)
    
    if(Q0.hat==0){
      tau <- 0
      if(sum(y>0)==1){
        warning("The Bootstrap community has only one species. Estimation is not robust.")
      }
    }else{
      tau <- U / t * (1 - Chat) / sum(y / t * (1 - y / t)^t)
    }
    
    pi <- y / t * (1 - tau * (1 - y / t)^t)
    pi.star <- c(pi, rep(U / t * (1 - Chat) / Q0, Q0))
    #   set.seed(456)
    y1 <- rbind(t, matrix(rbinom(length(pi.star) * B, t, pi.star), ncol = B))
    tmp <- which(colSums(y1)==t)
    if(length(tmp)>0) y1 <- y1[,-tmp]
    se <- sd(apply(y1, 2, function(y2) FunName(y2, datatype)))
    return(se)
  }
  
  if(datatype=="abundance"){
    BootstrapFun.abun(x=x, FunName, datatype, B)
  }else if(datatype=="incidence"){
    BootstrapFun.ince(y=x, FunName, datatype, B)
  }
}



#
#
###########################################
#' Estimation of Shannon entropy/diversity
#' 
#' \code{ChaoShannon}: estimation of Shannon entropy or transformed Shannon diversity based on the method proposed by Chao et al. (2013)
#' 
#' @param x a \code{matrix}, \code{data.frame} (species by sites), or \code{list} of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, 
#' then the first entry of the input data must be total number of sampling units, followed by species incidence frequencies. 
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param transform a \code{logical} constant to compute traditional Shannon entropy index (\code{transform=FALSE}) or the transformed Shannon diversity (\code{transform=TRUE}). 
#' @param conf a positive number \eqn{\le} 1 specifying the level of confidence interval. 
#' @param B an integer specifying the number of bootstrap replications.
#' @return A data.frame of observed Shannon entropy/diversity, estimate of entropy/diversity, s.e. and the associated confidence interval.
#' @seealso \code{\link{ChaoRichness}, \link{ChaoSimpson}} 
#' @examples 
#' data(spider)
#' ChaoShannon(spider$Girdled, datatype="abundance")
#' @references 
#' Chao, A., Wang, Y.T. & Jost, L. (2013) Entropy and the species accumulation curve: a novel entropy estimator via discovery rates of new species. Methods in Ecology and Evolution, 4, 1091-1100.
#'
#' @export


ChaoShannon <- function(x, datatype="abundance", transform=FALSE, conf=0.95, B=200) {
      
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid data type")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous data type")
  datatype <- match.arg(datatype, TYPE)
  
  if(datatype=="incidence_freq") datatype <- "incidence"
  
  if(datatype=="incidence_raw"){
    if(inherits(x, "list")){
      x <- lapply(x, as.incfreq)
    }else{
      x <- as.incfreq(x)
    }
    datatype <- "incidence"
  }
  
  if(!is.numeric(conf) || conf > 1 || conf < 0){
    warning("\"conf\"(confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!")
    conf <- 0.95
  }
  
  
  ChaoEntropyEstFun <- function(x, datatype){
    if(datatype == "abundance"){
      x <- x[x > 0]
      n <- sum(x)
      f1 <- sum(x==1)
      f2 <- sum(x==2)
      f0.hat <- ifelse(f2==0, (n-1)/n*f1*(f1-1)/2, (n-1)/n*f1^2/2/f2)  #estimation of unseen species via Chao1
      A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
      temp1 <- sum(x / n * (digamma(n) - digamma(x)))
      if(A == 1){
        temp2 <- 0
      }else{
        l <- 1:(n-1)
        temp2 <- f1 / n * (A)^(1-n) * (-log(1-A) - sum(1 / l * (A)^l))
      }
      temp2 <- ifelse(is.na(temp2), 0, temp2)
      est <- temp1 + temp2
    }else if(datatype == "incidence"){
      t <- x[1]
      x <- x[-1]
      x <- x[x > 0]
      U <- sum(x)
      Q1 <- sum(x==1)
      Q2 <- sum(x==2)
      Q0.hat <- ifelse(Q2==0, (t-1)/t*Q1*(Q1-1)/2, (t-1)/t*Q1^2/2/Q2)  #estimation of unseen species via Chao2
      A <- ifelse(Q1>0, t*Q0.hat/(t*Q0.hat+Q1), 1)
      temp1 <- sum(x / t * (digamma(t) - digamma(x)))
      
      if(A == 1){
        temp2 <- 0
      }else{
        r <- 1 : (t-1)
        temp2 <- Q1 / t * (A)^(1-t) * (-log(1-A) - sum(1 / r * (A)^r))
      }
      temp2 <- ifelse(is.na(temp2), 0, temp2)
      est <- t / U * (temp1 + temp2) + log(U / t)
    }
      return(est)
  }
  
  TransformEntropy <- function(x, datatype) exp(ChaoEntropyEstFun(x, datatype))
  
  ObsEntropy <- function(x, datatype){
    if(datatype=="abundance"){
      x <- x[x>0]
      p <- x/sum(x)
      -sum(p*log(p))
    }else if(datatype=="incidence"){
      #t <- x[1]
      y <- x[-1]
      y <- y[y>0]
      p <- y/sum(y)
      -sum(p*log(p))
    }
  } 
  
  z <- qnorm(1 - (1 - conf)/2)
  
  myFun <- function(x){
    if(!is.numeric(x) & !is.matrix(x) & !is.data.frame(x))
      stop("invalid data structure")
    if(is.matrix(x) | is.data.frame(x)){
      if (ncol(x) != 1 & nrow(x) != 1)
        stop("invalid data structure")
    }
    
    if(transform==TRUE){
      obs <- exp(ObsEntropy(x, datatype))
      est <- TransformEntropy(x, datatype)
      se <- BootstrapFun(x, TransformEntropy, datatype, B)
    }else{
      if(!is.logical(transform)){
        transform <- FALSE
        warning("transform must be a logical value, we use transform=FALSE to We use to calculate!")
      }
      obs <- ObsEntropy(x, datatype)
      est <- ChaoEntropyEstFun(x, datatype)
      se <- BootstrapFun(x, ChaoEntropyEstFun, datatype, B)
    }
    CI <- c(max(est - z * se, obs), est + z * se)
    out <- round(t(c(obs, est, se, CI)),3)
    out <- data.frame(out)
    colnames(out) <- c("Observed", "Estimator", "Est_s.e",
                     paste(conf*100, "% Lower", sep=""), paste(conf*100, "% Upper", sep=""))
    return(out)
  }
  
  if(inherits(x, "numeric")){
    out <- myFun(x)
  }else if(inherits(x, "integer")){
	out <- myFun(x)
  }else if(inherits(x, "list")){
    out <- do.call("rbind", lapply(x, myFun))
  } else if(inherits(x, "matrix") | inherits(x, "data.frame")){
    out <- do.call("rbind", apply(as.matrix(x), 2, myFun))  
  }
  return(out)
}

#
#
###########################################
#' Estimation of Gini-Simpson index or Simpson diversity
#' 
#' \code{ChaoSimpson}: estimation of Gini-Simpson index or the transformed Simpson diversity based on the methods proposed in Good (1953) and Chao et al. (2014)
#' 
#' @param x a \code{matrix}, \code{data.frame} (species by sites), or \code{list} of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, 
#' then the first entry of the input data must be total number of sampling units, followed by species incidence frequencies.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param transform a \code{logical} constant to compute traditional Gini-Simpson index (\code{transform=FALSE}) or the transformed Simpson diversity (\code{transform=TRUE}). 
#' @param conf a positive number \eqn{\le} 1 specifying the level of confidence interval. 
#' @param B an integer specifying the number of bootstrap replications.
#' @return a data.frame of observed Gini-Simpson index/diversity, index/diversity estimator, s.e. and the associated confidence interval.
#' @seealso \code{\link{ChaoRichness}, \link{ChaoShannon}}
#' @examples 
#' data(spider)
#' ChaoSimpson(spider$Girdled, datatype="abundance")
#' @references
#' Chao, A., Gotelli, N.J., Hsieh, T.C., Sander, E.L., Ma, K.H., Colwell, R.K. & Ellison, A.M. (2014) Rarefaction and extrapolation with Hill numbers: a framework for sampling and estimation in species diversity studies. Ecological Monographs, 84, 45-67.\cr\cr 
#' Good, I.J. (1953) The population frequencies of species and the estimation of population parameters. Biometrika, 40, 237-264.
#' 
#' @export

ChaoSimpson <- function(x, datatype="abundance", transform=FALSE, conf=0.95, B=200) {
    
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid data type")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous data type")
  datatype <- match.arg(datatype, TYPE)
  
  if(datatype=="incidence_freq") datatype <- "incidence"
  
  if(datatype=="incidence_raw"){
    if(inherits(x, "list")){
      x <- lapply(x, as.incfreq)
    }else{
      x <- as.incfreq(x)
    }
    datatype <- "incidence"
  }
  
  if(!is.numeric(conf) || conf > 1 || conf < 0){
    warning("\"conf\"(confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!")
    conf <- 0.95
  }
  Observed <- function(x, datatype){
    if(datatype=="abundance"){
      x <- x[x>0]
      p <- x/sum(x)
      1-sum(p^2)
    }else if(datatype=="incidence"){
      t <- x[1]
      y <- x[-1]
      y <- y[y>0]
      p <- y/t
      1-sum(p^2)/sum(p)^2
    }
  }
  MVUE <- function(x, datatype){
    if(datatype=="abundance"){
      x <- x[x>0]
      n <- sum(x)
      est <- 1-sum(choose(x,2)/choose(n,2))
      
    }else if(datatype=="incidence"){
      t <- x[1]
      y <- x[-1]
      y <- y[y>0]
      Q1 <- sum(y==1)
      a <- sum(choose(y,2)/choose(t,2))
      b <- (sum(y %*% t(y)) - sum(diag(y %*% t(y))))/t^2
      if(sum(y)!=Q1){
        est <- 1-a/(a+b)
      }else{
        est <- 1-2/t/(t+1)
      }
    }
    est
  }
  TransformSimpson <- function(x, datatype) {  
    if(max(x[-1]>1)){
      1/(1 - MVUE(x, datatype))
    }else {
      1/(1 - Observed(x, datatype))
    }
  }
  
  z <- qnorm(1 - (1 - conf)/2)
  
  myFun <- function(x){
    if(!is.numeric(x) & !is.matrix(x) & !is.data.frame(x))
      stop("invalid data structure")
    if(is.matrix(x) | is.data.frame(x)){
      if (ncol(x) != 1 & nrow(x) != 1)
        stop("invalid data structure")
    }
    if(transform==TRUE){
      obs <- 1/(1-Observed(x, datatype))
      est <- TransformSimpson(x, datatype)
      se <- BootstrapFun(x, TransformSimpson, datatype, B)
      CI <- c(max(est - z * se, obs), est + z * se)
    }else{
      if(!is.logical(transform)){
        transform <- FALSE
        warning("transform must be a logical value, we use transform=FALSE to We use to calculate!")
      }
      obs <- Observed(x, datatype)
      est <- MVUE(x, datatype)
      se <- BootstrapFun(x, MVUE, datatype, B)
      CI <- c(max(est - z * se, obs), min(est + z * se, 1))    
    }
    out <- round(t(c(obs, est, se, CI)),3)
    out <- data.frame(out)
    colnames(out) <- c("Observed", "Estimator", "Est_s.e.", paste(conf*100,"% Lower",sep=""), paste(conf*100,"% Upper",sep=""))  
    return(out)
  }
  
  if(inherits(x, "numeric")){
    out <- myFun(x)
  }else if(inherits(x, "integer")){
	out <- myFun(x)
  }else if(inherits(x, "list")){
    out <- do.call("rbind", lapply(x, myFun))
  } else if(inherits(x, "matrix") | inherits(x, "data.frame")){
    out <-do.call("rbind", apply(as.matrix(x), 2, myFun))  
  }
  return(out)
}

#' #' @export
#' #' @rdname ChaoRichness
#' ChaoSpecies <- ChaoRichness
#' 
#' #' @export
#' #' @rdname ChaoShannon
#' ChaoEntropy <- ChaoShannon
#' 
#' #' @export
#' #' @rdname ChaoSimpson
#' EstSimpson <- ChaoSimpson
