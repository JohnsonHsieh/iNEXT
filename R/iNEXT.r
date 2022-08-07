
#
#
###########################################
# Draw confidence band plot
# 
# \code{conf.reg} uses polygon to draw a confidence band plot
# 
# @param x a vector of estimator
# @param LCL a vector of lower bound
# @param UCL a vector of upwer bound
# @param ... further arguments to be passed to \code{polygon}
# @return a polygon plot
conf.reg=function(x,LCL,UCL, data=NULL,...) {
  x.sort <- order(x)
  x <- x[x.sort]
  LCL <- LCL[x.sort]
  UCL <- UCL[x.sort]
  polygon(c(x,rev(x)),c(LCL,rev(UCL)), ...)
}


#
#
###########################################
# Estimation of species reletive abundance distribution
# 
# \code{EstiBootComm.Ind} Estimation of species reletive abundance distribution to obtain bootstrap s.e.
# 
# @param Spec a vector of species abundances
# @return a vector of reltavie abundance
# @seealso \code{\link{EstiBootComm.Sam}}
# @examples 
# data(spider)
# EstiBootComm.Ind(spider$Girdled)
EstiBootComm.Ind <- function(Spec)
{
  Sobs <- sum(Spec > 0)   #observed species
  n <- sum(Spec)        #sample size
  f1 <- sum(Spec == 1)   #singleton 
  f2 <- sum(Spec == 2)   #doubleton
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  a <- f1/n*A
  b <- sum(Spec / n * (1 - Spec / n) ^ n)
  if(f0.hat==0){
    w <- 0
    if(sum(Spec>0)==1){
      warning("This site has only one species. Estimation is not robust.")
    }
  }else{
    w <- a / b      	#adjusted factor for rare species in the sample
  }
  Prob.hat <- Spec / n * (1 - w * (1 - Spec / n) ^ n)					#estimation of relative abundance of observed species in the sample
  Prob.hat.Unse <- rep(a/ceiling(f0.hat), ceiling(f0.hat))  	#estimation of relative abundance of unseen species in the sample
  return(sort(c(Prob.hat, Prob.hat.Unse), decreasing=TRUE))		  							#Output: a vector of estimated relative abundance
}


#
#
###########################################
# Estimation of species detection distribution
# 
# \code{EstiBootComm.Sam} Estimation of species detection distribution to obtain bootstrap s.e.
# 
# @param Spec a vector of species incidence, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
# @return a vector of estimated detection probability
# @seealso \code{\link{EstiBootComm.Sam}}
# @examples 
# data(ant)
# EstiBootComm.Sam(ant$h50m)
EstiBootComm.Sam <- function(Spec)
{
  nT <- Spec[1]
  Spec <- Spec[-1]
  Sobs <- sum(Spec > 0)   #observed species
  Q1 <- sum(Spec == 1) 	#singleton 
  Q2 <- sum(Spec == 2) 	#doubleton
  Q0.hat <- ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2)	#estimation of unseen species via Chao2
  A <- ifelse(Q1>0, nT*Q0.hat/(nT*Q0.hat+Q1), 1)
  a <- Q1/nT*A
  b <- sum(Spec / nT * (1 - Spec / nT) ^ nT)
  
  if(Q0.hat==0){
    w <- 0
    if(sum(Spec>0)==1){
      warning("This site has only one species. Estimation is not robust.")
    }
  }else{
    w <- a / b      	#adjusted factor for rare species in the sample
  }
  
  Prob.hat <- Spec / nT * (1 - w * (1 - Spec / nT) ^ nT)					#estimation of detection probability of observed species in the sample
  Prob.hat.Unse <- rep(a/ceiling(Q0.hat), ceiling(Q0.hat))  	#estimation of detection probability of unseen species in the sample
  return(sort(c(Prob.hat, Prob.hat.Unse), decreasing=TRUE))									#Output: a vector of estimated detection probability
}


#
#
###########################################
# iNterpolation and EXTrapolation of abundance-based Hill number
# 
# \code{TD.m.est} Estimation of interpolation and extrapolation of abundance-based Hill number with order q
# 
# @param x a vector of species abundances
# @param m a integer vector of rarefaction/extrapolation sample size
# @param qs a numerical vector of the order of Hill number
# @return a vector of estimated interpolation and extrapolation function of Hill number with order q
# @export
TD.m.est = function(x, m, qs){ ## here q is allowed to be a vector containing non-integer values.
  n <- sum(x)
  #xv_matrix = as.matrix(xv)
  ifi <- table(x);ifi <- cbind(i = as.numeric(names(ifi)),fi = ifi)
  obs <- Diversity_profile_MLE(x,qs)
  RFD_m <- RTD(ifi, n, n-1, qs)
  # RFD_m2 <- RTD(ifi, n, n-2, qs)
  # whincr <- which(RFD_m != RFD_m2)
  # Dn1 <- obs; Dn1[whincr] <- obs + (obs - RFD_m)^2/(RFD_m - RFD_m2)
  #asymptotic value
  asy <- Diversity_profile(x,qs)
  #beta
  beta <- rep(0,length(qs))
  # beta0plus <- which(asy != obs)
  # beta[beta0plus] <- (Dn1[beta0plus]-obs[beta0plus])/(asy[beta0plus]-obs[beta0plus])
  beta0plus <- which(asy != RFD_m)
  beta[beta0plus] <- (obs[beta0plus]-RFD_m[beta0plus])/(asy[beta0plus]-RFD_m[beta0plus])
  #Extrapolation, 
  ETD = function(m,qs){
    m = m-n
    out <- sapply(1:length(qs), function(i){
      if( qs[i] != 2) {
        obs[i]+(asy[i]-obs[i])*(1-(1-beta[i])^m)
      }else if( qs[i] == 2 ){
        1/ ((1/(n+m))+(1-1/(n+m))*sum(ifi[,2]*ifi[,1]/n*(ifi[,1]-1)/(n-1)) )
      } 
    })
    return(out)
  }
  Sub = function(m){
    if(m<n){
      if(m == round(m)) { mRTD[-1,mRTD[1,]==m] 
        } else { (ceiling(m)-m)*mRTD[-1,mRTD[1,]==floor(m)]+(m-floor(m))*mRTD[-1,mRTD[1,]==ceiling(m)] }
    }else if(m==n){
      obs
    }else{
      ETD(m,qs)
    }
  }
  
  if (sum(m < n) != 0) {
    int.m = sort(unique(c(floor(m[m<n]), ceiling(m[m<n]))))
    mRTD = rbind(int.m, sapply(int.m, function(k) RTD(ifi,n,k,qs)))
  }
  as.vector(t(sapply(m, Sub))) 
}

#
#
###########################################
# iNterpolation and EXTrapolation of incidence-based Hill number
# 
# \code{TD.m.est_inc} Estimation of interpolation and extrapolation of incidence-based Hill number
# 
# @param y a vector of species incidence-based frequency, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
# @param t_ a integer vector of rarefaction/extrapolation sample size
# @param qs a numerical vector of the order of Hill number
# @return a vector of estimated interpolation and extrapolation function of Hill number with order q
# @export
TD.m.est_inc <- function(y, t_, qs){
  nT <- y[1]
  y <- y[-1]
  y <- y[y > 0]
  U <- sum(y)
  #xv_matrix = as.matrix(xv)
  iQi <- table(y);iQi <- cbind(i = as.numeric(names(iQi)),Qi = iQi)
  obs <- Diversity_profile_MLE.inc(c(nT,y),qs)
  RFD_m <- RTD_inc(iQi, nT, nT-1, qs)
  # RFD_m2 <- RTD_inc(iQi, nT, nT-2, qs)
  # whincr <- which(RFD_m != RFD_m2)
  # Dn1 <- obs; Dn1[whincr] <- obs + (obs - RFD_m)^2/(RFD_m - RFD_m2)
  asy <- Diversity_profile.inc(c(nT,y),qs)
  beta <- rep(0,length(qs))
  # beta0plus <- which(asy != obs)
  # beta[beta0plus] <- (Dn1[beta0plus]-obs[beta0plus])/(asy[beta0plus]-obs[beta0plus])
  beta0plus <- which(asy != RFD_m)
  beta[beta0plus] <- (obs[beta0plus]-RFD_m[beta0plus])/(asy[beta0plus]-RFD_m[beta0plus])
  ETD = function(m,qs){
    m = m-nT
    out <- sapply(1:length(qs), function(i){
      if( qs[i] != 2) {
        obs[i]+(asy[i]-obs[i])*(1-(1-beta[i])^m)
      }else if( qs[i] == 2 ){
        1/ ((1/(nT+m))*(nT/U)+(1-1/(nT+m))*sum(iQi[,2]*iQi[,1]/(U^2)*(iQi[,1]-1)/(1-1/nT)) )
      } 
    })
    return(out)
  }
  Sub = function(m){
    if(m<nT){
      if(m == round(m)) { mRTD_inc[-1,mRTD_inc[1,]==m] 
      } else { (ceiling(m)-m)*mRTD_inc[-1,mRTD_inc[1,]==floor(m)]+(m-floor(m))*mRTD_inc[-1,mRTD_inc[1,]==ceiling(m)] }
    }else if(m==nT){
      obs
    }else{
      ETD(m,qs)
    }
  }
  
  if (sum(t_ < nT) != 0) {
    int.t_ = sort(unique(c(floor(t_[t_<nT]), ceiling(t_[t_<nT]))))
    mRTD_inc = rbind(int.t_, sapply(int.t_, function(k) RTD_inc(iQi,nT,k,qs)))
  }
  as.vector(t(sapply(t_, Sub)))
}

#
#
###############################################
# Abundance-based sample coverage
# 
# \code{Chat.Ind} Estimation of abundance-based sample coverage function
# 
# @param x a vector of species abundances
# @param m a integer vector of rarefaction/extrapolation sample size
# @return a vector of estimated sample coverage function
# @export
Chat.Ind <- function(x, m){
  x <- x[x>0]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  Sub <- function(m){
    #if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m < n) {
      xx <- x[(n-x)>=m]
      out <- 1-sum(xx / n * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }
    if(m == n) out <- 1-f1/n*A
    if(m > n) out <- 1-f1/n*A^(m-n+1)
    out
  }
  sapply(m, Sub)		
}


#
#
###############################################
# Incidence-based sample coverage
# 
# \code{Chat.Sam} Estimation of incidence-based sample coverage function
# 
# @param x a vector of species incidence-based frequency, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
# @param t a integer vector of rarefaction/extrapolation sample size
# @return a vector of estimated sample coverage function
# @export
Chat.Sam <- function(x, t){
  nT <- x[1]
  y <- x[-1]
  y <- y[y>0]
  U <- sum(y)
  Q1 <- sum(y == 1)
  Q2 <- sum(y == 2)
  Q0.hat <- ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2)  #estimation of unseen species via Chao2
  A <- ifelse(Q1>0, nT*Q0.hat/(nT*Q0.hat+Q1), 1)
  Sub <- function(t){
    if(t < nT) {
      yy <- y[(nT-y)>=t]
      out <- 1 - sum(yy / U * exp(lgamma(nT-yy+1)-lgamma(nT-yy-t+1)-lgamma(nT)+lgamma(nT-t)))     
    }
    #if(t < nT) out <- 1 - sum(y / U * exp(lchoose(nT - y, t) - lchoose(nT - 1, t)))
    if(t == nT) out <- 1 - Q1 / U * A
    if(t > nT) out <- 1 - Q1 / U * A^(t - nT + 1)
    out
  }
  sapply(t, Sub)		
}


#
#
###############################################
# iNterpolation and EXTrapolation of abundance-based Hill number
# 
# \code{iNEXT.Ind} Estimation of interpolation and extrapolation of abundance-based Hill number with order q
# 
# @param Spec a vector of species abundances
# @param q a numerical vector of the order of Hill number
# @param m a integer vector of rarefaction/extrapolation sample size, default is NULL. If m is not be specified, then the program will compute sample units due to endpoint and knots.
# @param endpoint a integer of sample size that is the endpoint for rarefaction/extrapolation. Default is double the original sample size.
# @param knots a number of knots of computation, default is 40
# @param se calculate bootstrap standard error and 95% confidence interval; default is TRUE
# @param nboot the number of bootstrap resampling times, default is 200
# @return a list of interpolation and extrapolation Hill number with specific order q (qD) and sample coverage (SC)
# @seealso \code{\link{iNEXT.Sam}}
# @examples
# data(spider)
# # q = 0 with specific endpoint
# iNEXT.Ind(spider$Girdled, q=0, endpoint=500)
# # q = 1 with specific sample size m and don't calculate standard error
# iNEXT.Ind(spider$Girdled, q=1, m=c(1, 10, 20, 50, 100, 200, 400, 600), se=FALSE)
iNEXT.Ind <- function(Spec, q=0, m=NULL, endpoint=2*sum(Spec), knots=40, se=TRUE, nboot=200, conf=0.95, unconditional_var = TRUE)
{
  qtile <- qnorm(1-(1-conf)/2)
  n <- sum(Spec)		  	#sample size
  if(is.null(m)) {
    if(endpoint <= n) {
      m <- floor(seq(1, endpoint, length=floor(knots)))
    } else {
      m <- c(floor(seq(1, sum(Spec)-1, length.out=floor(knots/2)-1)), sum(Spec), floor(seq(sum(Spec)+1, to=endpoint, length.out=floor(knots/2))))
    }
    m <- c(1, m[-1])
  } else if(is.null(m)==FALSE) {	
    if(max(m)>n | length(m[m==n])==0)  m <- c(m, n-1, n, n+1)
    m <- sort(m)
  }
  m <- unique(m)
  #====conditional on m====
  Dq.hat <- TD.m.est(Spec,m,q)
  C.hat <- Chat.Ind(Spec, m)
  #====unconditional====
  if(unconditional_var){
    goalSC <- unique(C.hat)
    Dq.hat_unc <- unique(invChat.Ind(x = Spec,q = q,C = goalSC))
    Dq.hat_unc$Method[round(Dq.hat_unc$m) == n] = "Observed"
  }

  if(se==TRUE & nboot > 1 & length(Spec) > 1) {
    Prob.hat <- EstiBootComm.Ind(Spec)
    Abun.Mat <- rmultinom(nboot, n, Prob.hat)
    
    ses_m <- apply(matrix(apply(Abun.Mat,2 ,function(x) TD.m.est(x, m, q)),
                        nrow = length(Dq.hat)),1,sd, na.rm=TRUE)
    
    ses_C_on_m <- apply(matrix(apply(Abun.Mat, 2, function(x) Chat.Ind(x, m)),nrow = length(m)),
                        1, sd, na.rm=TRUE)
    if(unconditional_var){
      ses_C <- apply(matrix(apply(Abun.Mat,2 ,function(x) invChat.Ind(x, q,goalSC)$qD),
                            nrow = nrow(Dq.hat_unc)),1,sd, na.rm=TRUE)
    }
  } else {
    ses_m <- rep(NA,length(Dq.hat))
    ses_C_on_m <- rep(NA,length(m))
    if(unconditional_var){
      ses_C <- rep(NA,nrow(Dq.hat_unc))
    }
  }
  out_m <- cbind("m"=rep(m,length(q)), "qD"=Dq.hat, "qD.LCL"=Dq.hat-qtile*ses_m,
                 "qD.UCL"=Dq.hat+qtile*ses_m,"SC"=rep(C.hat,length(q)), 
                 "SC.LCL"=C.hat-qtile*ses_C_on_m, "SC.UCL"=C.hat+qtile*ses_C_on_m)
  out_m <- data.frame(out_m)
  out_m$Method <- ifelse(out_m$m<n, "Rarefaction", ifelse(out_m$m==n, "Observed", "Extrapolation"))
  out_m$Order.q <- rep(q,each = length(m))
  id_m <- match(c("m", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL", "SC", "SC.LCL", "SC.UCL"), names(out_m), nomatch = 0)
  out_m <- out_m[, id_m]
  out_m$qD.LCL[out_m$qD.LCL<0] <- 0
  out_m$SC.LCL[out_m$SC.LCL<0] <- 0
  out_m$SC.UCL[out_m$SC.UCL>1] <- 1
  
  if(unconditional_var){
    out_C <- cbind(Dq.hat_unc,'qD.LCL' = Dq.hat_unc$qD-qtile*ses_C,
                   'qD.UCL' = Dq.hat_unc$qD+qtile*ses_C) 
    id_C <- match(c("SC","m", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL"), names(out_C), nomatch = 0)
    out_C <- out_C[, id_C]
    out_C$qD.LCL[out_C$qD.LCL<0] <- 0
  }else{
    out_C <- NULL
  }
  return(list(size_based = out_m, coverage_based = out_C))
}



#
#
###############################################
# iNterpolation and EXTrapolation of incidence-based Hill number
# 
# \code{iNEXT.Sam} Estimation of interpolation and extrapolation of incidence-based Hill number with order q
# 
# @param Spec a vector of species incidence-based frequency, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
# @param q a numerical vector of the order of Hill number
# @param t a integer vector of rarefaction/extrapolation sample size, default is NULL. If m is not be specified, then the program will compute sample units due to endpoint and knots.
# @param endpoint a integer of sample size that is the endpoint for rarefaction/extrapolation. Default is double the original sample size.
# @param knots a number of knots of computation, default is 40
# @param se calculate bootstrap standard error and 95% confidence interval; default is TRUE
# @param nboot the number of bootstrap resampling times, default is 200
# @return a list of interpolation and extrapolation Hill number with specific order q (qD) and sample coverage (SC)
# @seealso \code{\link{iNEXT.Ind}}
# @examples
# data(ant)
# # q = 0 with specific endpoint
# iNEXT.Sam(ant$h50m, q=0, endpoint=100)
# # q = 1 with specific sample size m and don't calculate standard error
# iNEXT.Sam(ant$h500m, q=1, t=round(seq(10, 500, length.out=20)), se=FALSE)
iNEXT.Sam <- function(Spec, t=NULL, q=0, endpoint=2*max(Spec), knots=40, se=TRUE, nboot=200, conf=0.95, unconditional_var = TRUE)
{
  qtile <- qnorm(1-(1-conf)/2)
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
    if(max(t)>nT | length(t[t==nT])==0)  t <- c(t, nT-1, nT, nT+1)
    t <- sort(t)
  }
  t <- unique(t)
  #====conditional on m====
  Dq.hat <- TD.m.est_inc(Spec,t,q)
  C.hat <- Chat.Sam(Spec, t)
  #====unconditional====
  # if(unconditional_var){
  #   goalSC <- unique(round(C.hat,4))
  #   goalSC[goalSC==1] <- 0.9999
  #   Dq.hat_unc <- unique(invChat.Sam(x = Spec,q = q,C = goalSC))
  # }
  if(unconditional_var){
    goalSC <- unique(C.hat)
    Dq.hat_unc <- unique(invChat.Sam(x = Spec,q = q,C = goalSC))
    Dq.hat_unc$Method[Dq.hat_unc$t == nT] = "Observed"
  }
  
  if(se==TRUE & nboot > 1 & length(Spec) > 2){
    Prob.hat <- EstiBootComm.Sam(Spec)
    Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
    Abun.Mat <- matrix(c(rbind(nT, Abun.Mat)),ncol=nboot)
    tmp <- which(colSums(Abun.Mat)==nT)
    if(length(tmp)>0) Abun.Mat <- Abun.Mat[,-tmp]
    if(ncol(Abun.Mat)==0){
      out <- cbind("t"=t, "qD"=Dq.hat, "SC"=C.hat)
      warning("Insufficient data to compute bootstrap s.e.")
    }else{		
      ses_m <- apply(matrix(apply(Abun.Mat,2 ,function(y) TD.m.est_inc(y, t, q)),
                            nrow = length(Dq.hat)),1,sd, na.rm=TRUE)
      
      ses_C_on_m <- apply(matrix(apply(Abun.Mat, 2, function(y) Chat.Sam(y, t)),nrow = length(t)),
                          1, sd, na.rm=TRUE)
      if(unconditional_var){
        ses_C <- apply(matrix(apply(Abun.Mat,2 ,function(y) invChat.Sam(y, q,goalSC)$qD),
                              nrow = nrow(Dq.hat_unc)),1,sd, na.rm=TRUE)
      }
    }
  }else {
    ses_m <- rep(NA,length(Dq.hat))
    ses_C_on_m <- rep(NA,length(t))
    if(unconditional_var){
      ses_C <- rep(NA,nrow(Dq.hat_unc))
    }
  }
  
  out_m <- cbind("t"=rep(t,length(q)), "qD"=Dq.hat, "qD.LCL"=Dq.hat-qtile*ses_m,
                 "qD.UCL"=Dq.hat+qtile*ses_m,"SC"=rep(C.hat,length(q)), 
                 "SC.LCL"=C.hat-qtile*ses_C_on_m, "SC.UCL"=C.hat+qtile*ses_C_on_m)
  out_m <- data.frame(out_m)
  out_m$Method <- ifelse(out_m$t<nT, "Rarefaction", ifelse(out_m$t==nT, "Observed", "Extrapolation"))
  out_m$Order.q <- rep(q,each = length(t))
  id_m <- match(c("t", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL", "SC", "SC.LCL", "SC.UCL"), names(out_m), nomatch = 0)
  out_m <- out_m[, id_m]
  out_m$qD.LCL[out_m$qD.LCL<0] <- 0
  out_m$SC.LCL[out_m$SC.LCL<0] <- 0
  out_m$SC.UCL[out_m$SC.UCL>1] <- 1
  
  if(unconditional_var){
    out_C <- cbind(Dq.hat_unc,'qD.LCL' = Dq.hat_unc$qD-qtile*ses_C,
                   'qD.UCL' = Dq.hat_unc$qD+qtile*ses_C) 
    id_C <- match(c("SC","t", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL"), names(out_C), nomatch = 0)
    out_C <- out_C[, id_C]
    out_C$qD.LCL[out_C$qD.LCL<0] <- 0
  }else{
    out_C <- NULL
  }
  return(list(size_based = out_m, coverage_based = out_C))
}


#
#
###############################################
#' iNterpolation and EXTrapolation of Hill numbers
#' 
#' \code{iNEXT}: Interpolation and extrapolation of Hill number with order q
#' 
#' @param x a \code{matrix}, \code{data.frame} (species by sites), or \code{list} of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units in each column or list. 
#' @param q a number or vector specifying the diversity order(s) of Hill numbers.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
# @param rowsum a logical variable to check if the input object is raw data (species by sites matrix, \code{rowsum=FALSE}) or iNEXT default input (abundance counts or incidence frequencies, \code{rowsum=TRUE}).
#' @param size an integer vector of sample sizes (number of individuals or sampling units) for which diversity estimates will be computed. 
#' If NULL, then diversity estimates will be computed for those sample sizes determined by the specified/default \code{endpoint} and \code{knots} .
#' @param endpoint an integer specifying the sample size that is the \code{endpoint} for rarefaction/extrapolation. 
#' If NULL, then \code{endpoint} \code{=} double the reference sample size.
#' @param knots an integer specifying the number of equally-spaced \code{knots} (say K, default is 40) between size 1 and the \code{endpoint};
#' each knot represents a particular sample size for which diversity estimate will be calculated.  
#' If the \code{endpoint} is smaller than the reference sample size, then \code{iNEXT()} computes only the rarefaction esimates for approximately K evenly spaced \code{knots}. 
#' If the \code{endpoint} is larger than the reference sample size, then \code{iNEXT()} computes rarefaction estimates for approximately K/2 evenly spaced \code{knots} between sample size 1 and the reference sample size, and computes extrapolation estimates for approximately K/2 evenly spaced \code{knots} between the reference sample size and the \code{endpoint}.
#' @param se a logical variable to calculate the bootstrap standard error and \code{conf} confidence interval.
#' @param conf a positive number < 1 specifying the level of confidence interval; default is 0.95.
#' @param nboot an integer specifying the number of replications; default is 50.
#' @importFrom reshape2 dcast
#' @importFrom stats rmultinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom stats rbinom
#' @importFrom stats optimize
#' @importFrom stats quantile
#' @return a list of three objects: \code{$DataInfo} for summarizing data information; 
#' \code{$iNextEst} for showing diversity estimates for rarefied and extrapolated samples along with related statistics;
#' and \code{$AsyEst} for showing asymptotic diversity estimates along with related statistics.\cr 
#' 
#' NOTE: From version 3.0.0, \code{$iNextEst} has been expanded to include \code{$size_based} and \cr 
#'       \code{$coverage_based} to provide two types of confidence intervals. 
#'
#' @examples
#' ## example for abundance based data (list of vector)
#' data(spider)
#' out1 <- iNEXT(spider, q=c(0,1,2), datatype="abundance")
#' out1$DataInfo # showing basic data information.
#' out1$AsyEst # showing asymptotic diversity estimates.
#' out1$iNextEst$size_based 
#' # showing diversity estimates with rarefied and extrapolated samples; 
#' # confidence limits are obtained for fixed sample size.
#' 
#' out1$iNextEst$coverage_based 
#' # showing diversity estimates with rarefied and extrapolated samples;
#' # confidence limits are obtained for fixed sample coverage.
#'
#' ## example for abundance based data (data.frame)
#' data(bird)
#' out2 <- iNEXT(bird, q=0, datatype="abundance")
#' out2
#' 
#' \dontrun{
#' ## example for incidence frequencies based data (list of data.frame)
#' data(ant)
#' t <- round(seq(10, 500, length.out=20))
#' out3 <- iNEXT(ant$h500m, q=1, datatype="incidence_freq", size=t, se=FALSE)
#' out3$iNextEst
#' }
#' @export
#' 
iNEXT <- function(x, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)
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
  
  Fun <- function(x, q, assem_name){
    x <- as.numeric(unlist(x))
    unconditional_var <- TRUE
    if(datatype == "abundance"){
      if(sum(x)==0) stop("Zero abundance counts in one or more sample sites")
      out <- iNEXT.Ind(Spec=x, q=q, m=size, endpoint=ifelse(is.null(endpoint), 2*sum(x), endpoint), knots=knots, se=se, nboot=nboot, conf=conf,unconditional_var)
    }
    if(datatype == "incidence"){
      t <- x[1]
      y <- x[-1]
      if(t>sum(y)){
        warning("Insufficient data to provide reliable estimators and associated s.e.") 
      }
      if(sum(x)==0) stop("Zero incidence frequencies in one or more sample sites")
      
      out <- iNEXT.Sam(Spec=x, q=q, t=size, endpoint=ifelse(is.null(endpoint), 2*max(x), endpoint), knots=knots, se=se, nboot=nboot, conf=conf)  
    }
    if(unconditional_var){
      out <- lapply(out, function(out_) cbind(Assemblage = assem_name, out_))
    }else{
      out[[1]] <- cbind(Assemblage = assem_name, out[[1]])
    }
    out
  }
  
  if(!inherits(q, "numeric"))
    stop("invlid class of order q, q should be a postive value/vector of numeric object")
  if(min(q) < 0){
    warning("ambigous of order q, we only compute postive q")
    q <- q[q >= 0]
  }
  
  z <- qnorm(1-(1-0.95)/2)
  if(class_x=="numeric" | class_x=="integer" | class_x=="double"){
    out <- Fun(x,q,'Assemblage1')
    
    out <- list(size_based = out[[1]],
                coverage_based = out[[2]])
    
    index <- AsyD(x = x,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq')
                     ,nboot = 100,conf = 0.95)
    LCL <- index$qD.LCL[index$method=='Estimated']
    UCL <- index$qD.UCL[index$method=='Estimated']
    index <- dcast(index,formula = Order.q~method,value.var = 'qD')
    index <- cbind(index[,-1],se = (UCL - index$Estimated)/z,LCL,UCL)
    index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
    colnames(index) <- c("Observed","Estimator","Est_s.e.","95% Lower","95% Upper")
    # index <- rbind(as.matrix(ChaoSpecies(x, datatype, conf)), 
    #                as.matrix(ChaoEntropy(x, datatype, transform=TRUE, conf)),
    #                as.matrix(EstSimpson(x, datatype, transform=TRUE, conf)))
    rownames(index) <- c("Species Richness", "Shannon diversity", "Simpson diversity")
  
    
  }else if(class_x=="matrix" | class_x=="data.frame"){
    if(is.null(colnames(x))){
      colnames(x) <- sapply(1:ncol(x), function(i) paste0('assemblage',i))
    }
    out <- lapply(1:ncol(x), function(i) {
      tmp <- Fun(x[,i],q,colnames(x)[i])
      tmp
    })
    out <- list(size_based = do.call(rbind,lapply(out,  function(out_){out_[[1]]})),
                coverage_based = do.call(rbind,lapply(out,  function(out_){out_[[2]]})))
    
    index <- AsyD(x = x,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq'),nboot = 100,conf = 0.95)
    index = index[order(index$Assemblage),]
    LCL <- index$qD.LCL[index$method=='Estimated']
    UCL <- index$qD.UCL[index$method=='Estimated']
    index <- dcast(index,formula = Assemblage+Order.q~method,value.var = 'qD')
    index <- cbind(index,se = (UCL - index$Estimated)/z,LCL,UCL)
    index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
    index$Order.q <- c('Species richness','Shannon diversity','Simpson diversity')
    # arr <- array(0, dim = c(3, 5, ncol(x)))
    # arr[1,,] <- t(as.matrix(ChaoSpecies(x, datatype, conf)))
    # arr[2,,] <- t(as.matrix(ChaoEntropy(x, datatype, transform=TRUE, conf)))
    # arr[3,,] <- t(as.matrix(EstSimpson(x, datatype, transform=TRUE, conf)))  
    # dimnames(arr)[[3]] <- names(x)
    # dimnames(arr)[[1]] <- c("Species richness", "Shannon diversity", "Simpson diversity")
    # dimnames(arr)[[2]] <- c("Observed", "Estimator", "Est_s.e.", "Lower_CI", "Upper_CI")
    # index <- ftable(arr, row.vars = c(3,1))
    # index <- dcast(as.data.frame(index), formula = Var1+Var2~Var3, value.var = "Freq")
    colnames(index) <- c("Assemblage", "Diversity", "Observed", "Estimator", "s.e.", "LCL", "UCL")
    
    
  }else if(class_x=="list"){
    if(is.null(names(x))){
      names(x) <- sapply(1:length(x), function(i) paste0('assemblage',i))
    }
    out <- lapply(1:length(x), function(i) {
      tmp <- Fun(x[[i]],q,names(x)[i])
      tmp
    })
    
    out <- list(size_based = do.call(rbind,lapply(out,  function(out_){out_[[1]]})),
                coverage_based = do.call(rbind,lapply(out,  function(out_){out_[[2]]})))
    
    index <- AsyD(x = x,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq'),nboot = 100,conf = 0.95)
    index = index[order(index$Assemblage),]
    LCL <- index$qD.LCL[index$method=='Estimated']
    UCL <- index$qD.UCL[index$method=='Estimated']
    index <- dcast(index,formula = Assemblage+Order.q~method,value.var = 'qD')
    index <- cbind(index,se = (UCL - index$Estimated)/z,LCL,UCL)
    index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
    index$Order.q <- c('Species richness','Shannon diversity','Simpson diversity')
    # arr <- array(0, dim = c(3, 5, length(x)))
    # arr[1,,] <- t(as.matrix(ChaoSpecies(x, datatype, conf)))
    # arr[2,,] <- t(as.matrix(ChaoEntropy(x, datatype, transform=TRUE, conf)))
    # arr[3,,] <- t(as.matrix(EstSimpson(x, datatype, transform=TRUE, conf)))  
    # dimnames(arr)[[3]] <- names(x)
    # dimnames(arr)[[1]] <- c("Species richness", "Shannon diversity", "Simpson diversity")
    # dimnames(arr)[[2]] <- c("Observed", "Estimator", "Est_s.e.", "Lower_CI", "Upper_CI")
    # index <- ftable(arr, row.vars = c(3,1))
    # index <- dcast(as.data.frame(index), formula = Var1+Var2~Var3, value.var = "Freq")
    colnames(index) <- c("Assemblage", "Diversity", "Observed", "Estimator", "s.e.", "LCL", "UCL")
  }else{
    stop("invalid class of x, x should be a object of numeric, matrix, data.frame, or list")
  }
  out$size_based$Assemblage <- as.character(out$size_based$Assemblage)
  out$coverage_based$Assemblage <- as.character(out$coverage_based$Assemblage)
  info <- DataInfo(x, datatype)
  
  
  z <- list("DataInfo"=info, "iNextEst"=out, "AsyEst"=index)
  class(z) <- c("iNEXT")
  return(z)
  }


#
#
###########################################
# Estimation of the rank of species relative abundance distribution or detection probability
# 
# \code{EstDis} Estimation of the rank of species relative abundance distribution or detection probability to obtain bootstrap s.e.
# 
# @param x a vector of species abundances or incidence frequencies. If \code{datatype = "incidence"}, then the input format of first entry should be total number of sampling units, and followed by species incidence frequency.
# @param datatype the data type of input data. That is individual-based abundance data (\code{datatype = "abundance"}) or sample-based incidence data (\code{datatype = "incidence"}).
# @return a vector of the rank of estimated relative abundance distribution or detection probability
# @examples 
# data(spider)
# EstDis(spider$Girdled, datatype="abundance")
# data(ant)
# EstDis(ant$h50m, datatype="incidence_freq")
# @export
EstDis <- function(x, datatype=c("abundance", "incidence")){
  datatype <- match.arg(datatype)
  if(datatype == "abundance") out <- EstiBootComm.Ind(Spec=x)                                                      
  if(datatype == "incidence") out <- EstiBootComm.Sam(Spec=x)
  out 
}

#
#
###############################################
# Asymptotic diversity q profile 
AsyD <- function(x, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95){
  if(datatype == "incidence"){
    stop('datatype="incidence" was no longer supported, 
         please try datatype = "incidence_freq".')  
  }
  if(datatype == "incidence_raw"){
    stop('datatype="incidence_raw" was no longer supported, 
         please try datatype = "incidence_freq".')  
  }
  TYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  if(!inherits(q, "numeric"))
    stop("invlid class of order q, q should be a postive value/vector of numeric object")
  if(min(q) < 0){
    warning("ambigous of order q, we only compute postive q")
    q <- q[q >= 0]
  }
  if(nboot < 0 | round(nboot) - nboot != 0)
    stop("Please enter non-negative integer for nboot.")
  if(conf < 0 | conf > 1)
    stop("Please enter value between zero and one for confident interval.")
  
  if (inherits(x, "data.frame") | inherits(x, "matrix")){
    datalist <- lapply(1:ncol(x), function(i) x[,i])
    if(is.null(colnames(x))) names(datalist) <-  paste0("data",1:ncol(x)) else names(datalist) <- colnames(x)
    x <- datalist
  } else if (inherits(x, "numeric") | inherits(x, "integer") | inherits(x, "double")) {
    x <- list(data = x)
  }
  
  if(datatype=="abundance"){
    out <- lapply(1:length(x),function(i){
      dq <- c(Diversity_profile(x[[i]],q),Diversity_profile_MLE(x[[i]],q))
      if(nboot > 1){
        Prob.hat <- EstiBootComm.Ind(x[[i]])
        Abun.Mat <- rmultinom(nboot, sum(x[[i]]), Prob.hat)
        
        error <- qnorm(1-(1-conf)/2) * 
          apply(apply(Abun.Mat, 2, function(xb) c(Diversity_profile(xb, q),Diversity_profile_MLE(xb,q))), 1, sd, na.rm=TRUE)
        
      } else {error = NA}
      out <- data.frame("Order.q" = rep(q,2), "qD" = dq,"qD.LCL" = dq - error, "qD.UCL" = dq + error,
                        "Assemblage" = names(x)[i], "method" = rep(c("Estimated","Empirical"),each = length(q)))
      out$qD.LCL[out$qD.LCL<0] <- 0
      out
    })
    out <- do.call(rbind,out)
  }else if(datatype=="incidence_freq"){
    out <- lapply(1:length(x),function(i){
      dq <- c(Diversity_profile.inc(x[[i]],q),Diversity_profile_MLE.inc(x[[i]],q))
      if(nboot > 1){
        nT <- x[[i]][1]
        Prob.hat <- EstiBootComm.Sam(x[[i]])
        Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
        Abun.Mat <- matrix(c(rbind(nT, Abun.Mat)),ncol=nboot)
        tmp <- which(colSums(Abun.Mat)==nT)
        if(length(tmp)>0) Abun.Mat <- Abun.Mat[,-tmp]
        if(ncol(Abun.Mat)==0){
          error = 0
          warning("Insufficient data to compute bootstrap s.e.")
        }else{		
          error <- qnorm(1-(1-conf)/2) * 
            apply(apply(Abun.Mat, 2, function(yb) c(Diversity_profile.inc(yb, q),Diversity_profile_MLE.inc(yb,q))), 1, sd, na.rm=TRUE)
        }
      } else {error = NA}
      out <- data.frame("Order.q" = rep(q,2), "qD" = dq,"qD.LCL" = dq - error, "qD.UCL" = dq + error,
                        "Assemblage" = names(x)[i],"method" = rep(c("Estimated","Empirical"),each = length(q)))
      out$qD.LCL[out$qD.LCL<0] <- 0
      out
    })
    out <- do.call(rbind,out)
  }
  return(out)
}


##
##
###########################################
## Example individual-based data, spiders abundance data collected by Sackett et al. (2011)
##
##
# Girdled <- c(46, 22, 17, 15, 15, 9, 8, 6, 6, 4, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
# Logged <- c(88, 22, 16, 15, 13, 10, 8, 8, 7, 7, 7, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
# spider <- list(Girdled=Girdled, Logged=Logged)

##
##
###########################################
## Example sample-based data, tropical ant species data collected by Longino and Colwell (2011)
## Note that first cell is number of total samples, and others are species incidence-based frequency.
##

## 50m
# y50 <- c(599,rep(1,49),rep(2,23),rep(3,18),rep(4,14),rep(5,9),rep(6,10),rep(7,4),
#          rep(8,8),rep(9,6),rep(10,2),rep(11,1),12,12,13,13,rep(14,5),15,15,
#          rep(16,4),17,17,17,18,18,19,19,20,20,20,21,22,23,23,25,27,27,29,30,30,
#          31,33,39,40,43,46,46,47,48,51,52,52,56,56,58,58,61,61,65,69,72,77,79,82,
#          83,84,86,91,95,97,98,98,106,113,124,126,127,128,129,129,182,183,186,195,
#          222,236,263,330)

##500m
# y500 <- c(230,rep(1,71),rep(2,34),rep(3,12),rep(4,14),rep(5,9),rep(6,11),rep(7,8),
#           rep(8,4),rep(9,7),rep(10,5),rep(11,2),12,12,12,13,13,13,13,14,14,15,
#           16,16,17,17,17,17,18,19,20,21,21,23,24,25,25,25,26,27,30,31,31,32,32,
#           33,34,36,37,38,38,38,38,39,39,41,42,43,44,45,46,47,49,52,52,53,54,56,
#           60,60,65,73,78,123,131,133)

##1070m
# y1070 <- c(150,rep(1,28),rep(2,16),rep(3,13),rep(4,3),rep(5,1),rep(6,3),rep(7,6),
#            rep(8,1),rep(9,1),rep(10,1),rep(11,4),12,12,12,13,13,13,13,14,15,
#            16,16,16,16,18,19,19,21,22,23,24,25,25,25,26,30,31,31,31,32,34,36,
#            38,39,43,43,45,45,46,54,60,68,74,80,96,99)
##1500m
# y1500 <- c(200,rep(1,13),rep(2,4),rep(3,2),rep(4,2),rep(5,4),rep(6,2),rep(9,4),
#            rep(11,2),rep(17,2),18,19,23,23,24,25,25,25,29,30,32,33,43,50,53,
#            73,74,76,79,113,144)

##2000m
# y2000=c(200,1,2,2,3,4,8,8,13,15,19,23,34,59,80)
# 
# ant <- list(h50m=y50, h500m=y500, h1070m=y1070, h1500m=y1500, h2000m=y2000)


#' @useDynLib iNEXT, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL