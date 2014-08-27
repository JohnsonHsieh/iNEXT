
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
conf.reg=function(x,LCL,UCL,...) {
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
  w <- a / b    		#adjusted factor for rare species in the sample
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
  w <- a / b  			#adjusted factor for rare species in the sample
  Prob.hat <- Spec / nT * (1 - w * (1 - Spec / nT) ^ nT)					#estimation of detection probability of observed species in the sample
  Prob.hat.Unse <- rep(a/ceiling(Q0.hat), ceiling(Q0.hat))  	#estimation of detection probability of unseen species in the sample
  return(sort(c(Prob.hat, Prob.hat.Unse), decreasing=TRUE))									#Output: a vector of estimated detection probability
}


#
#
###########################################
# iNterpolation and EXTrapolation of abundance-based Hill number
# 
# \code{Dqhat.Ind} Estimation of interpolation and extrapolation of abundance-based Hill number with order q
# 
# @param x a vector of species abundances
# @param q a numerical value of the order of Hill number
# @param m a integer vector of rarefaction/extrapolation sample size
# @return a vector of estimated interpolation and extrapolation function of Hill number with order q
# @export
Dqhat.Ind <- function(x, q, m){
	x <- x[x > 0]
	n <- sum(x)
	
	fk.hat <- function(x, m){
		x <- x[x > 0]
		n <- sum(x)
		if(m <= n){
			Sub <- function(k)	sum(exp(lchoose(x, k) + lchoose(n - x, m -k) - lchoose(n, m)))
			sapply(1:m, Sub)
		}
		
		else {
			f1 <- sum(x == 1)
			f2 <- sum(x == 2)
			A <- ifelse(f2 > 0, (n-1)*f1/((n-1)*f1+2*f2), (n-1)*f1/((n-1)*f1+2))
			C.hat <- 1 - f1 / n * A
			p.hat <- x / n * C.hat			
			Sub <- function(k)	sum((choose(m, k) * p.hat^k * (1 - p.hat)^(m - k)) / (1 - (1 - p.hat)^n))
			sapply(1:m, Sub)
		}
	}
	
	D0.hat <- function(x, m){
		x <- x[x > 0]
		n <- sum(x)
		Sub <- function(m){
			if(m <= n){
				Fun <- function(x){
					if(x <= (n - m)) exp(lgamma(n - x + 1) + lgamma(n - m + 1) - lgamma(n - x - m + 1) - lgamma(n + 1))
					else 0
				}
				sum(1 - sapply(x, Fun))
			}
			else {
				Sobs <- sum(x > 0)
				f1 <- sum(x == 1)
				f2 <- sum(x == 2)
				f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)	#estimation of unseen species via Chao1
				A <- n*f0.hat/(n*f0.hat+f1)
				ifelse(f1 ==0, Sobs ,Sobs + f0.hat * (1 - A ^ (m - n)))	
			}
		}
		sapply(m, Sub)
	}
	
	D1.hat <- function(x, m){
		x <- x[x > 0]
		n <- sum(x)
		Sub <- function(m){
		  if(m < n){
			k <- 1:m
			exp(-sum(k / m * log(k / m) * fk.hat(x, m)))
		  }
		  else{
			#UE=sum(sapply(1:(n-1),function(k){sum(1/k*x/n*exp(lchoose(n-x,k)-lchoose(n-1,k)))}))
			UE <- sum(x/n*(digamma(n)-digamma(x)))
			f1 <- sum(x == 1)
			f2 <- sum(x == 2)
			A <- 1 - ifelse(f2 > 0, (n-1)*f1/((n-1)*f1+2*f2), (n-1)*f1/((n-1)*f1+2))
			#A=2*sum(x==2)/((n-1)*sum(x==1)+2*sum(x==2))
			B=sum(x==1)/n*(1-A)^(-n+1)*(-log(A)-sum(sapply(1:(n-1),function(k){1/k*(1-A)^k})))
			H.hat <- UE+B
			Hn.hat <- -sum(x / n * log(x / n))
			w <- (m - n) / m
			
			exp(w * H.hat + (1 - w) * Hn.hat)
			
		  }
		}
		sapply(m, Sub)
	}
	
	D2.hat <- function(x, m){
		Sub <- function(m) 1 / (1 / m + (1 - 1 / m) * sum(x * (x - 1) / n / (n - 1)))
		sapply(m, Sub)
	}
	
	Dq.hat <- function(x, m){
		Sub <- function(m){
			k <- 1:m
			sum( (k / m)^q * fk.hat(x, m))^(1 / (1 - q))
		}
		sapply(m, Sub)
	}
	if(q == 0) D0.hat(x, m)
	else if(q == 1) D1.hat(x, m)
	else if(q == 2) D2.hat(x, m)
	else Dq.hat(x, m)
}


#
#
###########################################
# iNterpolation and EXTrapolation of incidence-based Hill number
# 
# \code{Dqhat.Sam} Estimation of interpolation and extrapolation of incidence-based Hill number
# 
# @param y a vector of species incidence-based frequency, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
# @param q a numerical value of the order of Hill number
# @param t a integer vector of rarefaction/extrapolation sample size
# @return a vector of estimated interpolation and extrapolation function of Hill number with order q
# @export
Dqhat.Sam <- function(y, q, t){

	nT <- y[1]
	y <- y[-1]
	y <- y[y > 0]
	U <- sum(y)
	
	Qk.hat <- function(y, nT, t){
		if(t <= nT){
			Sub <- function(k)	sum(exp(lchoose(y, k) + lchoose(nT - y, t - k) - lchoose(nT, t)))
			sapply(1:t, Sub)
		}
		
		else {
			p.hat <- EstiBootComm.Sam(c(nT, y))
			Sub <- function(k)	sum((choose(t, k) * p.hat^k * (1 - p.hat)^(t - k)) / (1 - (1 - p.hat)^T))
			sapply(1:t, Sub)
		}
	}
	
	D0.hat <- function(y, nT, t){
		Sub <- function(t){
			if(t <= nT){
				Fun <- function(y){
					if(y <= (nT - t)) exp(lgamma(nT - y + 1) + lgamma(nT - t + 1) - lgamma(nT - y - t + 1) - lgamma(nT + 1))
					else 0
				}
				sum(1 - sapply(y, Fun))
			}
			else {
				Sobs <- sum(y > 0)
				Q1 <- sum(y==1)
				Q2 <- sum(y==2)
				Q0.hat <- ifelse(Q2 == 0,  (nT-1)/nT* Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2)	#estimation of unseen species via Chao2
				A <- nT*Q0.hat/(nT*Q0.hat+Q1)
				ifelse(Q1 ==0, Sobs ,Sobs + Q0.hat * (1 - A ^ (t - nT)))	
			}
		}
		sapply(t, Sub)
	}
	
	D1.hat <- function(y, nT, t){
		U <- sum(y)
		Sub <- function(t){
		  if(t < nT){
			k <- 1:t	
			Ut.hat <- t / nT * U
			exp(-sum(k / Ut.hat * log(k / Ut.hat) * Qk.hat(y, nT, t)))
		  }
		  else {
			UE <- sum(y / nT * (digamma(nT) - digamma(y)))
			Q1 <- sum(y == 1)
			Q2 <- sum(y == 2)
			A <- 1 - ifelse(Q2 > 0, (nT-1)*Q1/((nT-1)*Q1+2*Q2), (nT-1)*Q1/((nT-1)*Q1+2))
			B <- sum(y==1)/nT*(1-A)^(-nT+1)*(-log(A)-sum(sapply(1:(nT-1),function(k){1/k*(1-A)^k})))
			H.hat <- UE+B
			H.hat <- nT/U*H.hat-log(nT/U)
		  
			Hn.hat <- -sum(y / U * log(y / U))
			w <- (t - nT) / t
			exp(w * H.hat + (1 - w) * Hn.hat)
		  }
		}
		sapply(t, Sub)
	}
	
	D2.hat <- function(y, nT, t){
		U <- sum(y)
		Sub <- function(t) 1 / (1 / t * nT / U + (1 - 1 / t) * sum(y * (y - 1) / U^2 / (1 - 1 / nT)))
		sapply(t, Sub)
	}
	
	Dq.hat <- function(y, nT, t){
		U <- sum(y)
		Sub <- function(t){
			k <- 1:t
			Ut.hat <- U * t / nT
			sum( (k / Ut.hat)^q * Qk.hat(y, nT, t))^(1 / (1 - q))
		}
		sapply(t, Sub)
	}
	if(q == 0) D0.hat(y, nT, t)
	else if(q == 1) D1.hat(y, nT, t)
	else if(q == 2) D2.hat(y, nT, t)
	else Dq.hat(y, nT, t)
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
		if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
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
		if(t < nT) out <- 1 - sum(y / U * exp(lchoose(nT - y, t) - lchoose(nT - 1, t)))
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
# @param q a numeric value, the order of Hill number 
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
iNEXT.Ind <- function(Spec, q=0, m=NULL, endpoint=2*sum(Spec), knots=40, se=TRUE, nboot=200)
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
		if(max(m)>n & length(m[m==n])==0)  m <- c(m, n, n+1)
		m <- sort(m)
	}
	
	Dq.hat <- Dqhat.Ind(Spec, q, m)
	C.hat <- Chat.Ind(Spec, m)
	
	if(se==TRUE & nboot > 0 & length(Spec) > 1) {
		Prob.hat <- EstiBootComm.Ind(Spec)
		Abun.Mat <- rmultinom(nboot, n, Prob.hat)
	
		error <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(x) Dqhat.Ind(x, q, m)), 1, sd, na.rm=TRUE)
		left  <- Dq.hat - error
		right <- Dq.hat + error
	
		error.C <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(x) Chat.Ind(x, m)), 1, sd, na.rm=TRUE)
		left.C  <- C.hat - error.C
		right.C <- C.hat + error.C
		out <- cbind("m"=m, "qD"=Dq.hat, "qD.95.LCL"=left, "qD.95.UCL"=right, "SC"=C.hat, "SC.95.LCL"=left.C, "SC.95.UCL"=right.C)
	} else {
		out <- cbind("m"=m, "qD"=Dq.hat, "SC"=C.hat)
	}
	out <- data.frame(out)
  out$method <- ifelse(out$m<n, "interpolated", ifelse(out$m==n, "observed", "extrapolated"))
  out$order <- q
	id <- match(c("m", "method", "order", "qD", "qD.95.LCL", "qD.95.UCL", "SC", "SC.95.LCL", "SC.95.UCL"), names(out), nomatch = 0)
  out <- out[, id]
	return(out)
}



#
#
###############################################
# iNterpolation and EXTrapolation of incidence-based Hill number
# 
# \code{iNEXT.Sam} Estimation of interpolation and extrapolation of incidence-based Hill number with order q
# 
# @param Spec a vector of species incidence-based frequency, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
# @param q a numeric value, the order of Hill number 
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
iNEXT.Sam <- function(Spec, t=NULL, q=0, endpoint=2*max(Spec), knots=40, se=TRUE, nboot=200)
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
		if(max(t)>nT & length(t[t==nT])==0)  t <- c(t, nT, nT+1)
		t <- sort(t)
	}
	
	Dq.hat <- Dqhat.Sam(Spec, q, t)
	C.hat <- Chat.Sam(Spec, t)
	
	if(se==TRUE & nboot > 0 & length(Spec) > 2){
		Prob.hat <- EstiBootComm.Sam(Spec)
		Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
		Abun.Mat <- matrix(c(rbind(nT, Abun.Mat)),ncol=nboot)
		
		error <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(y) Dqhat.Sam(y, q, t)), 1, sd, na.rm=TRUE)
		left  <- Dq.hat - error
		right <- Dq.hat + error
		left[left<=0] <- 0
				
		error.C <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(y) Chat.Sam(y, t)), 1, sd, na.rm=TRUE)
		left.C  <- C.hat - error.C
		right.C <- C.hat + error.C
		left.C[left.C<=0] <- 0
    right.C[right.C>=1] <- 1
    
		out <- cbind("t"=t, "qD"=Dq.hat, "qD.95.LCL"=left, "qD.95.UCL"=right, "SC"=C.hat, "SC.95.LCL"=left.C, "SC.95.UCL"=right.C)
	} else {
		out <- cbind("t"=t, "qD"=Dq.hat, "SC"=C.hat)
	}
	out <- data.frame(out)
	out$method <- ifelse(out$t<nT, "interpolated", ifelse(out$t==nT, "observed", "extrapolated"))
	out$order <- q
  id <- match(c("t", "method", "order", "qD", "qD.95.LCL", "qD.95.UCL", "SC", "SC.95.LCL", "SC.95.UCL"), names(out), nomatch = 0)
	out <- out[, id]
	return(out)
}


#
#
###############################################
#' iNterpolation and EXTrapolation of Hill number
#' 
#' \code{iNEXT} Estimation of interpolation and extrapolation of Hill number with order q
#' 
#' @param x a vector of species abundance or incidence frequency. If \code{datatype = "incidence"}, then the input format of first entry should be total number of sampling units, and followed by species incidence frequency. 
#' Note that in iNEXT version 2.0, the input arguments allow the class of \code{\link{numeric}}, \code{\link{matrix}}, \code{\link{data.frame}} (species by sites), or \code{\link{list}}.
#' @param q a numeric value, the order of Hill number .
#' @param datatype the data type of input data. That is individual-based abundance data (\code{datatype = "abundance"}) or sample-based incidence data (\code{datatype = "incidence"}).
#' @param size an integer vector of rarefaction/extrapolation sample size (number of individuals or sampling units), default is NULL. If \code{size} is not be specified, it would compute rarefaction/extrapolation by endpoint and knots.
#' @param endpoint a integer of sample size that is the endpoint for rarefaction/extrapolation. Default is double reference sample size.
#' @param knots a number of knots (say K, default is 40) specifying separate sample size between 1 and \code{endpoint}. 
#' If \code{endpoint} is smaller than reference sample size, then \code{iNEXT()} compute rarefaction part only and divided by approximately equall spaced K. 
#' If \code{endpoint} os larger than reference sample size, then \code{iNEXT()} will compute approximately K/2 equally spaced for rarefaction part and divided extrapolation part as approximately K/2 equally spaced between reference sample size and \code{endpoint}.
#' @param se calculate bootstrap standard error and show 95\% confidence interval; default is TRUE.
#' @param nboot the number of bootstrap resampling times, default is 50.
#' @return a list of interpolation and extrapolation Hill number with specific order q (qD) and sample coverage (SC)
#' @examples
#' data(spider)
#' z <- iNEXT(spider, q=0, datatype="abundance")
#' 
#' data(ant)
#' iNEXT(ant$h500m, q=1, datatype="incidence", size=round(seq(10, 500, length.out=20)), se=FALSE)
#' @export
#' 
iNEXT <- function(x, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, nboot=50)
{
  TYPE <- c("abundance", "incidence")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  
  Fun <- function(x, q){
    if(datatype == "abundance")
      out <- iNEXT.Ind(Spec=x, q=q, m=size, endpoint=ifelse(is.null(endpoint), 2*sum(x), endpoint), knots=knots, se=se, nboot=nboot)
    if(datatype == "incidence")
      out <- iNEXT.Sam(Spec=x, q=q, t=size, endpoint=ifelse(is.null(endpoint), 2*max(x), endpoint), knots=knots, se=se, nboot=nboot)  
    out
  }
  
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a postive value/vector of numeric object")
  if(min(q) < 0){
    warning("ambigous of order q, we only compute postive q")
    q <- q[q >= 0]
  }
  
  if(class(x)=="numeric"){
    out <- do.call("rbind", lapply(q, function(q) Fun(x, q)))
    out[,-(1:3)] <- round(out[,-(1:3)],3)
    index <- rbind(as.matrix(ChaoSpecies(x, datatype)), 
                   as.matrix(ChaoEntropy(x, datatype, transform=TRUE)),
                   as.matrix(EstSimpson(x, datatype, transform=TRUE)))
    rownames(index) <- c("Species Richness", "Exponential Entropy", "Inverse Simpson")
    
  }else if(class(x)=="matrix" | class(x)=="data.frame"){
  
    out <- apply(as.matrix(x), 2, function(x){
      tmp <- do.call("rbind", lapply(q, function(q) Fun(x,q)))
      tmp[,-(1:3)] <- round(tmp[,-(1:3)],3)
      tmp
      })
    arr <- array(0, dim = c(3, 5, ncol(x)))
    arr[1,,] <- t(as.matrix(ChaoSpecies(x, datatype)))
    arr[2,,] <- t(as.matrix(ChaoEntropy(x, datatype, transform=TRUE)))
    arr[3,,] <- t(as.matrix(EstSimpson(x, datatype, transform=TRUE)))  
    dimnames(arr)[[3]] <- names(x)
    dimnames(arr)[[1]] <- c("Species Richness", "Exponential Entropy", "Inverse Simpson")
    dimnames(arr)[[2]] <- c("Observed", "Estimator", "Est_s.e.", "95% Lower", "95% Upper")
    index <- ftable(arr, row.vars = c(3,1))
    
  }else if(class(x)=="list"){
    out <- lapply(x, function(x) {
      tmp <- do.call("rbind", lapply(q, function(q) Fun(x,q)))
      tmp[,-(1:3)] <- round(tmp[,-(1:3)],3)
      tmp
      })
    
    arr <- array(0, dim = c(3, 5, length(x)))
    arr[1,,] <- t(as.matrix(ChaoSpecies(x, datatype)))
    arr[2,,] <- t(as.matrix(ChaoEntropy(x, datatype, transform=TRUE)))
    arr[3,,] <- t(as.matrix(EstSimpson(x, datatype, transform=TRUE)))  
    dimnames(arr)[[3]] <- names(x)
    dimnames(arr)[[1]] <- c("Species Richness", "Exponential Entropy", "Inverse Simpson")
    dimnames(arr)[[2]] <- c("Observed", "Estimator", "Est_s.e.", "95% Lower", "95% Upper")
    index <- ftable(arr, row.vars = c(3,1))
  }else{
    stop("invlid class of x, x should be a object of numeric, matrix, data.frame, or list")
  }
  
  info <- DataInfo(x, datatype)

  
  z <- list("DataInfo"=info, "BasicIndex"=index,"Accumulation"=out)
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
# @param x a vector of species abundance or incidence frequency. If \code{datatype = "incidence"}, then the input format of first entry should be total number of sampling units, and followed by species incidence frequency.
# @param datatype the data type of input data. That is individual-based abundance data (\code{datatype = "abundance"}) or sample-based incidence data (\code{datatype = "incidence"}).
# @return a vector of the rank of estimated relative abundance distribution or detection probability
# @examples 
# data(spider)
# EstDis(spider$Girdled, datatype="incidence")
# data(ant)
# EstDis(ant$h50m, datatype="incidence")
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
#' ggplot2 extension for an iNEXT Object
#' 
#' \code{ggiNEXT} the \code{\link{ggplot2}} extension for \code{\link{iNEXT}} Object
#' @param x a \code{iNEXT} object computed by \code{\link{iNEXT}}
#' @param type three different plotting \code{type = c(1, 2, 3)}; \code{1} means to plot number of individuals or number of samples to diversity;
#'              \code{2} means to plot number of individuals or number of samples to sample coverage; and \code{3} means to plot sample coverage to diversity.                 
#' @param se display confidence interval around estimated accumulation curve
#' @param facet.var display subsets of the dataset in different panels with four choices: \code{facet.var = c("none", "order", "site", "both")} where \code{"none"} means do not split any lay out panels in a grid; \code{"order"} means split lay out panels by different orders q; \code{"site"} means split lay out panels by different sites; and \code{"both"} means split lay out panels by order and sites.              
#' @param color.var display subsets of the dataset in different colors with four choices: \code{color.var = c("none", "order", "site", "both")} where \code{"none"} means do not split any color; \code{"order"} means split colors by different orders q; \code{"site"} means split colors by different sites; and \code{"both"} means split colors by order and sites.  
#' @return a ggplot object
#' @examples
#' data(spider)
#' # single abundance-based data
#' out1 <- iNEXT(spider$Girdled, q=0, datatype="abundance")
#' ggiNEXT(x=out1, type=1)
#' ggiNEXT(x=out1, type=2)
#' ggiNEXT(x=out1, type=3)
#' 
#' # single incidence-based data with multiple order q
#' data(ant)
#' y <- iNEXT(ant$h500m, q=c(0,1,2), datatype="incidence", size=round(seq(10, 500, length.out=20)), se=FALSE)
#' ggiNEXT(y, se=FALSE)
#' 
#' # multiple abundance-based data with multiple order q
#' z <- iNEXT(spider, q=c(0,1,2), datatype="abundance")
#' ggiNEXT(z, facet.var="site", color.var="order")
#' ggiNEXT(z, facet.var="both", color.var="both")
#' @export

ggiNEXT <- function(x, type=1, se=TRUE, facet.var="none", color.var="order"){
  if(class(x) != "iNEXT") 
    stop("invalid object class")
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

    
  z <- x$Accumulation
  if(class(z) == "list"){
    z <- data.frame(do.call("rbind", z) ,site=rep(names(z), sapply(z, nrow)))
    rownames(z) <- NULL
  }
  
  if("qD.95.LCL" %in% names(z) == FALSE & se) {
    warning("invalid se setting, the iNEXT object do not consist confidence interval")
    se <- FALSE
  }else if("qD.95.LCL" %in% names(z) & se) {
    se <- TRUE
  }else{
    se <- FALSE
  }
  
  
  if(type==1L) {
    z$x <- z[,1]
    z$y <- z$qD
    if(se){
      z$y.lwr <- z$qD.95.LCL
      z$y.upr <- z$qD.95.UCL
    }
  }else if(type==2L){
    z$x <- z[,1]
    z$y <- z$SC
    if(se){
      z$y.lwr <- z$SC.95.LCL
      z$y.upr <- z$SC.95.UCL
    }
  }else if(type==3L){
    z$x <- z$SC
    z$y <- z$qD
    if(se){
      z$y.lwr <- z$qD.95.LCL
      z$y.upr <- z$qD.95.UCL
    }
  }
  
  y <- method <- site <- y.lwr <- y.upr <- NULL
  if(color.var=="none"){
    if(levels(factor(z$order))>1 & "site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT object consists multiple sites and orders, change setting as both")
      color.var <- "both"
      z$col <- paste(z$site, z$order, sep="-")
    }else if("site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT object consists multiple orders, change setting as order")
      color.var <- "site"
      z$col <- z$site
    }else if(levels(factor(z$order))>1){
      warning("invalid color.var setting, the iNEXT object consists multiple sites, change setting as site")
      color.var <- "order"
      z$col <- z$order
    }else{
      z$col <- rep(1, nrow(z))
    }
  }else if(color.var=="order"){     
    z$col <- z$order
  }else if(color.var=="site"){
    if(!"site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- order
    }
    z$col <- z$site
  }else if(color.var=="both"){
    if(!"site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- order
    }
    z$col <- paste(z$site, z$order, sep="-")
  }


  
  if("site"%in%names(z)){
    g <- ggplot(z, aes(x=x, y=y, colour=factor(col))) + geom_point(aes(x=x, y=y, colour=factor(col)), size=5, data=subset(z, method=="observed"))
  }else{
    g <- ggplot(z, aes(x=x, y=y, colour=factor(col), shape=site)) + geom_point(aes(x=x, y=y, colour=factor(col), shape=site), size=5, data=subset(z, method=="observed"))
  }
  
  g <- g + geom_line(aes(linetype=factor(method, c("interpolated", "extrapolated"), c("interpolation", "extrapoltion"))), size=1.5) +        
    labs(x="Number of sampling units", y="Species diversity") +
    guides(linetype=guide_legend(title="Method"), 
           colour=guide_legend(title="Order"), 
           fill=guide_legend(title="Order"), 
           shape=guide_legend(title="Site")) + 
    theme(legend.position = "bottom", text=element_text(size=18)) 
  
  if(se)
    g <- g + geom_ribbon(aes(ymin=y.lwr, ymax=y.upr, fill=factor(col), colour=NULL), alpha=0.2)
  
  
  if(facet.var=="order"){
    if(length(levels(factor(z$order))) == 1){
      warning("invalid facet.var setting, the iNEXT object do not consist multiple orders.")      
    }else{
      g <- g + facet_grid(.~ order)
      if(color.var=="both"){
        g <- g + guides(colour=guide_legend(title="Site-Order", ncol=length(levels(factor(z$order))), byrow=TRUE),
                        fill=guide_legend(title="Site-Order"))
      }
    }
  }
  
  if(facet.var=="site"){
    if(!"site"%in%names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites.")
    }else{
      g <- g + facet_grid(.~site)
      if(color.var=="both"){
        g <- g + guides(colour=guide_legend(title="Site-Order", nrow=length(levels(factor(z$order)))),
                        fill=guide_legend(title="Site-Order"))
      }
    }
  }
  
  if(facet.var=="both"){
    if(length(levels(factor(z$order))) == 1 | !"site"%in%names(z)){
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites or orders.")
    }else{
      g <- g + facet_wrap(site~order) 
      if(color.var=="both"){
        g <- g +  guides(colour=guide_legend(title="Site-Order", nrow=length(levels(factor(z$site))), byrow=TRUE),
                         fill=guide_legend(title="Site-Order"))
      }
    }
  }
  return(g)
  
}



##
##
###########################################
## Example individual-based data, spiders abundance data collected by Sackett et al. (2011)
##
##
Girdled <- c(46, 22, 17, 15, 15, 9, 8, 6, 6, 4, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
Logged <- c(88, 22, 16, 15, 13, 10, 8, 8, 7, 7, 7, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
spider <- list(Girdled=Girdled, Logged=Logged)

##
##
###########################################
## Example sample-based data, tropical ant species data collected by Longino and Colwell (2011)
## Note that first cell is number of total samples, and others are species incidence-based frequency.
##

## 50m
y50 <- c(599,rep(1,49),rep(2,23),rep(3,18),rep(4,14),rep(5,9),rep(6,10),rep(7,4),
		rep(8,8),rep(9,6),rep(10,2),rep(11,1),12,12,13,13,rep(14,5),15,15,
		rep(16,4),17,17,17,18,18,19,19,20,20,20,21,22,23,23,25,27,27,29,30,30,
		31,33,39,40,43,46,46,47,48,51,52,52,56,56,58,58,61,61,65,69,72,77,79,82,
		83,84,86,91,95,97,98,98,106,113,124,126,127,128,129,129,182,183,186,195,
		222,236,263,330)
	
##500m
y500 <- c(230,rep(1,71),rep(2,34),rep(3,12),rep(4,14),rep(5,9),rep(6,11),rep(7,8),
		rep(8,4),rep(9,7),rep(10,5),rep(11,2),12,12,12,13,13,13,13,14,14,15,
		16,16,17,17,17,17,18,19,20,21,21,23,24,25,25,25,26,27,30,31,31,32,32,
		33,34,36,37,38,38,38,38,39,39,41,42,43,44,45,46,47,49,52,52,53,54,56,
		60,60,65,73,78,123,131,133)

##1070m
y1070 <- c(150,rep(1,28),rep(2,16),rep(3,13),rep(4,3),rep(5,1),rep(6,3),rep(7,6),
		rep(8,1),rep(9,1),rep(10,1),rep(11,4),12,12,12,13,13,13,13,14,15,
		16,16,16,16,18,19,19,21,22,23,24,25,25,25,26,30,31,31,31,32,34,36,
		38,39,43,43,45,45,46,54,60,68,74,80,96,99)
##1500m
y1500 <- c(200,rep(1,13),rep(2,4),rep(3,2),rep(4,2),rep(5,4),rep(6,2),rep(9,4),
		rep(11,2),rep(17,2),18,19,23,23,24,25,25,25,29,30,32,33,43,50,53,
		73,74,76,79,113,144)

##2000m
y2000=c(200,1,2,2,3,4,8,8,13,15,19,23,34,59,80)

ant <- list(h50m=y50, h500m=y500, h1070m=y1070, h1500m=y1500, h2000m=y2000)
