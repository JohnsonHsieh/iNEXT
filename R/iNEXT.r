#
#
###########################################
#' Draw confidence band plot
#' 
#' \code{conf.reg} uses polygon to draw a confidence band plot
#' 
#' @param x a vector of estimator
#' @param LCL a vector of lower bound
#' @param UCL a vector of upwer bound
#' @param ... further arguments to be passed to \code{polygon}
#' @return a polygon plot
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
#' Estimation of species reletive abundance distribution
#' 
#' \code{EstiBootComm.Ind} Estimation of species reletive abundance distribution to obtain bootstrap s.e.
#' 
#' @param Spec a vector of species abundances
#' @return a vector of reltavie abundance
#' @seealso \code{\link{EstiBootComm.Sam}}
#' @examples 
#' data(spider)
#' EstiBootComm.Ind(spider$Girdled)
#' @export
EstiBootComm.Ind <- function(Spec)
{
  Sobs <- sum(Spec > 0)   #observed species
  n <- sum(Spec)        #sample size
  f1 <- sum(Spec == 1)   #singleton 
  f2 <- sum(Spec == 2)   #doubleton
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- n*f0.hat/(n*f0.hat+f1)
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
#' Estimation of species detection distribution
#' 
#' \code{EstiBootComm.Sam} Estimation of species detection distribution to obtain bootstrap s.e.
#' 
#' @param Spec a vector of species incidence, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
#' @return a vector of estimated detection probability
#' @seealso \code{\link{EstiBootComm.Sam}}
#' @examples 
#' data(ant)
#' EstiBootComm.Sam(ant$h50m)
#' @export
EstiBootComm.Sam <- function(Spec)
{
  nT <- Spec[1]
  Spec <- Spec[-1]
  Sobs <- sum(Spec > 0)   #observed species
  Q1 <- sum(Spec == 1) 	#singleton 
  Q2 <- sum(Spec == 2) 	#doubleton
  Q0.hat <- ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2)	#estimation of unseen species via Chao2
  A <- nT*Q0.hat/(nT*Q0.hat+Q1)
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
#' iNterpolation and EXTrapolation of abundance-based Hill number
#' 
#' \code{Dqhat.Ind} Estimation of interpolation and extrapolation of abundance-based Hill number with order q
#' 
#' @param x a vector of species abundances
#' @param q a numerical value of the order of Hill number
#' @param m a integer vector of rarefaction/extrapolation sample size
#' @return a vector of estimated interpolation and extrapolation function of Hill number with order q

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
#' iNterpolation and EXTrapolation of incidence-based Hill number
#' 
#' \code{Dqhat.Sam} Estimation of interpolation and extrapolation of incidence-based Hill number
#' 
#' @param y a vector of species incidence-based frequency, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
#' @param q a numerical value of the order of Hill number
#' @param t a integer vector of rarefaction/extrapolation sample size
#' @return a vector of estimated interpolation and extrapolation function of Hill number with order q
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
		  if(t < T){
			k <- 1:t	
			Ut.hat <- t / nT * U
			exp(-sum(k / Ut.hat * log(k / Ut.hat) * Qk.hat(y, nT, t)))
		  }
		  else {
			UE <- sum(y / nT * (digamma(nT) - digamma(y)))
			Q1 <- sum(y == 1)
			Q2 <- sum(y == 2)
			A <- 1 - ifelse(Q2 > 0, (nT-1)*Q1/((nT-1)*Q1+2*Q2), (nT-1)*Q1/((nT-1)*Q1+2))
			B=sum(y==1)/nT*(1-A)^(-nT+1)*(-log(A)-sum(sapply(1:(nT-1),function(k){1/k*(1-A)^k})))
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
#' Abundance-based sample coverage
#' 
#' \code{Chat.Ind} Estimation of abundance-based sample coverage function
#' 
#' @param x a vector of species abundances
#' @param m a integer vector of rarefaction/extrapolation sample size
#' @return a vector of estimated sample coverage function
Chat.Ind <- function(x, m){
	x <- x[x>0]
	n <- sum(x)
	f1 <- sum(x == 1)
	f2 <- sum(x == 2)
	f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
	A <- n*f0.hat/(n*f0.hat+f1)
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
#' Incidence-based sample coverage
#' 
#' \code{Chat.Sam} Estimation of incidence-based sample coverage function
#' 
#' @param x a vector of species incidence-based frequency, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
#' @param t a integer vector of rarefaction/extrapolation sample size
#' @return a vector of estimated sample coverage function
Chat.Sam <- function(x, t){
	nT <- x[1]
	y <- x[-1]
	y <- y[y>0]
	U <- sum(y)
	Q1 <- sum(y == 1)
	Q2 <- sum(y == 2)
	Q0.hat <- ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2)  #estimation of unseen species via Chao2
	A <- nT*Q0.hat/(nT*Q0.hat+Q1)
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
#' iNterpolation and EXTrapolation of abundance-based Hill number
#' 
#' \code{iNEXT.Ind} Estimation of interpolation and extrapolation of abundance-based Hill number with order q
#' 
#' @param Spec a vector of species abundances
#' @param q a numeric value, the order of Hill number 
#' @param m a integer vector of rarefaction/extrapolation sample size, default is NULL. If m is not be specified, then the program will compute sample units due to endpoint and knots.
#' @param endpoint a integer of sample size that is the endpoint for rarefaction/extrapolation. Default is double the original sample size.
#' @param Knots a number of knots of computation, default is 40
#' @param se calculate bootstrap standard error and 95% confidence interval; default is TRUE
#' @param nboot the number of bootstrap resampling times, default is 200
#' @return a list of interpolation and extrapolation Hill number with specific order q (qD) and sample coverage (SC)
#' @seealso \code{\link{iNEXT.Sam}}
#' @examples
#' data(spider)
#' # q = 0 with specific endpoint
#' iNEXT.Ind(spider$Girdled, q=0, endpoint=500)
#' # q = 1 with specific sample size m and don't calculate standard error
#' iNEXT.Ind(spider$Girdled, q=1, m=c(1, 10, 20, 50, 100, 200, 400, 600), se=FALSE)
#' @export
iNEXT.Ind <- function(Spec, q=0, m=NULL, endpoint=2*sum(Spec), Knots=40, se=TRUE, nboot=200)
{
	n <- sum(Spec)		  	#sample size
	if(is.null(m)) {
		if(endpoint <= n) {
			m <- floor(seq(1, endpoint, length=floor(Knots)))
		} else {
			m <- c(floor(seq(1, sum(Spec)-1, length.out=floor(Knots/2)-1)), sum(Spec), floor(seq(sum(Spec)+1, to=endpoint, length.out=floor(Knots/2))))
		}
		m <- c(1, m[-1])
	} else if(is.null(m)==FALSE) {	
		if(max(m)>n & length(m[m==n])==0)  m <- c(m, n)
		m <- sort(m)
	}
	
	Dq.hat <- Dqhat.Ind(Spec, q, m)
	C.hat <- Chat.Ind(Spec, m)
	
	if(se==TRUE & nboot > 0) {
		Prob.hat <- EstiBootComm.Ind(Spec)
		Abun.Mat <- rmultinom(nboot, n, Prob.hat)
	
		error <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(x) Dqhat.Ind(x, q, m)), 1, sd, na.rm=TRUE)
		left  <- Dq.hat - error
		right <- Dq.hat + error
	
		error.C <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(x) Chat.Ind(x, m)), 1, sd, na.rm=TRUE)
		left.C  <- C.hat - error.C
		right.C <- C.hat + error.C
		out <- cbind("m"=m, "qD"=Dq.hat, "qD.95%LCL"=left, "qD.95%UCL"=right, "SC"=C.hat, "SC.95%LCL"=left.C, "SC.95%UCL"=right.C)
	} else {
		out <- cbind("m"=m, "qD"=Dq.hat, "SC"=C.hat)
	}
	out <- data.frame(out)
	if(max(m) > n){
		out.int <- out[m<=n,]
		out.ext <- out[m>n,]
		return(list("interpolation"=out.int, "extrapolation"=out.ext))
	} else {
		return(list("interpolation"=out))
	}
}



#
#
###############################################
#' iNterpolation and EXTrapolation of abundance-based Hill number
#' 
#' \code{iNEXT.Sam} Estimation of interpolation and extrapolation of incidence-based Hill number with order q
#' 
#' @param Spec a vector of species incidence-based frequency, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
#' @param q a numeric value, the order of Hill number 
#' @param t a integer vector of rarefaction/extrapolation sample size, default is NULL. If m is not be specified, then the program will compute sample units due to endpoint and knots.
#' @param endpoint a integer of sample size that is the endpoint for rarefaction/extrapolation. Default is double the original sample size.
#' @param Knots a number of knots of computation, default is 40
#' @param se calculate bootstrap standard error and 95% confidence interval; default is TRUE
#' @param nboot the number of bootstrap resampling times, default is 200
#' @return a list of interpolation and extrapolation Hill number with specific order q (qD) and sample coverage (SC)
#' @seealso \code{\link{iNEXT.Ind}}
#' @examples
#' data(ant)
#' # q = 0 with specific endpoint
#' iNEXT.Sam(ant$h50m, q=0, endpoint=100)
#' # q = 1 with specific sample size m and don't calculate standard error
#' iNEXT.Sam(ant$h500m, q=1, t=round(seq(10, 500, length.out=20)), se=FALSE)
#' @export
iNEXT.Sam <- function(Spec, t=NULL, q=0, endpoint=2*max(Spec), Knots=40, se=TRUE, nboot=200)
{
	if(which.max(Spec)!=1) {return(warning("Wrong input format!"))}
	nT <- Spec[1]
	if(is.null(t)) {
		if(endpoint <= nT) {
			t <- floor(seq(1, endpoint, length.out=floor(Knots)))
		} else {
			t <- c(floor(seq(1, nT-1, length.out=floor(Knots/2)-1)), nT, floor(seq(nT+1, to=endpoint, length.out=floor(Knots/2))))
		}
		t <- c(1, t[-1])
	} else if(is.null(t)==FALSE) {	
		if(max(t)>nT & length(t[t==nT])==0)  t <- c(t, nT)
		t <- sort(t)
	}
	
	Dq.hat <- Dqhat.Sam(Spec, q, t)
	C.hat <- Chat.Sam(Spec, t)
	
	if(se==TRUE & nboot > 0){
		Prob.hat <- EstiBootComm.Sam(Spec)
		Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
		Abun.Mat <- matrix(c(rbind(nT, Abun.Mat)),ncol=nboot)
		
		error <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(y) Dqhat.Sam(y, q, t)), 1, sd, na.rm=TRUE)
		left  <- Dq.hat - error
		right <- Dq.hat + error
				
		error.C <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(y) Chat.Sam(y, t)), 1, sd, na.rm=TRUE)
		left.C  <- C.hat - error.C
		right.C <- C.hat + error.C
		out <- cbind("t"=t, "qD"=Dq.hat, "qD.95%LCL"=left, "qD.95%UCL"=right, "SC"=C.hat, "SC.95%LCL"=left.C, "SC.95%UCL"=right.C)
	} else {
		out <- cbind("t"=t, "qD"=Dq.hat, "SC"=C.hat)
	}
	out <- data.frame(out)
	if(max(t) > nT){
		out.int <- out[t<=nT,]
		out.ext <- out[t>nT,]
		return(list("interpolation"=out.int, "extrapolation"=out.ext))
	} else {
		return(list("interpolation"=out))
	}
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
