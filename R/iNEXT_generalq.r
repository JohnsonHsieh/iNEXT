#Asymptotic diversity
Diversity_profile <- function(x,q){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  sortx = sort(unique(x))
  tab = table(x)
  Sub_q012 <- function(q){
    if(q==0){
      length(x) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }else if(q==1){
      A <- sum(tab*sortx/n*(digamma(n)-digamma(sortx)))
      B <- D1_2nd(n,f1,p1)
      exp(A+B)
    }else if(abs(q-round(q))==0){
      A <- sum(tab[sortx>=q]*exp(lchoose(sortx[sortx>=q],q)-lchoose(n,q)))
      A^(1/(1-q))
    }
  }
  ans <- rep(0,length(q))
  q_part1 = which(abs(q-round(q))==0)
  if(length(q_part1)>0){
    ans[q_part1] <- sapply(q[q_part1], Sub_q012)
  }
  q_part2 <- which(!abs(q-round(q))==0)
  if(length(q_part2)>0){
    ans[q_part2] <- Dq(ifi = cbind(i = sortx, fi = tab),n = n,qs = q[q_part2],f1 = f1,A = p1)
  }
  ans
}
Diversity_profile.inc <- function(data,q){
  T = data[1]
  Yi = data[-1]
  Yi <- Yi[Yi!=0]
  U <- sum(Yi)
  Q1 <- sum(Yi==1)
  Q2 <- sum(Yi==2)
  Sobs <- length(Yi)
  A <- AA.inc(data)
  Q0hat <- ifelse(Q2 == 0, (T - 1) / T * Q1 * (Q1 - 1) / 2, (T - 1) / T * Q1 ^ 2/ 2 / Q2)
  B <- sapply(q,function(q) ifelse(A==1,0,(Q1/T)*(1-A)^(-T+1)*(A^(q-1)-sum(sapply(c(0:(T-1)),function(r) choose(q-1,r)*(A-1)^r)))))
  qD <- (U/T)^(q/(q-1))*(qDFUN(q,Yi,T) + B)^(1/(1-q))
  qD[which(q==0)] = Sobs+Q0hat
  yi <- Yi[Yi>=1 & Yi<=(T-1)]
  delta <- function(i){
    (yi[i]/T)*sum(1/c(yi[i]:(T-1)))
  }
  if(sum(q %in% 1)>0){
    C_ <- ifelse(A==1,0,(Q1/T)*(1-A)^(-T+1)*(-log(A)-sum(sapply(c(1:(T-1)),function(r) (1-A)^r/r))))
    qD[which(q==1)] <- exp((T/U)*( sum(sapply(c(1:length(yi)),function(i) delta(i))) + C_)+log(U/T))
  }
  return(qD)
}
#diversity at reference sample size
Diversity_profile_MLE <- function(x,q){
  p <- x[x>0]/sum(x)
  Sub <- function(q){
    if(q==0) sum(p>0)
    else if(q==1) exp(-sum(p*log(p)))
    else exp(1/(1-q)*log(sum(p^q)))
  }
  sapply(q, Sub)
}
Diversity_profile_MLE.inc <- function(data,q){
  Yi = data[-1]
  U = sum(Yi)
  Yi <- Yi[Yi!=0]
  ai <- Yi/U
  qD = qD_MLE(q,ai)
  qD[which(q==1)] <- exp(-sum(ai*log(ai)))
  return(qD)
}
#iNEXT
TD.m.est = function(x, m, qs){ ## here q is allowed to be a vector containing non-integer values.
  n <- sum(x)
  #xv_matrix = as.matrix(xv)
  ifi <- table(x);ifi <- cbind(i = as.numeric(names(ifi)),fi = ifi)
  obs <- Diversity_profile_MLE(x,qs)
  RFD_m <- RTD(ifi, n, n-1, qs)
  RFD_m2 <- RTD(ifi, n, n-2, qs)
  whincr <- which(RFD_m != RFD_m2)
  Dn1 <- obs; Dn1[whincr] <- obs + (obs - RFD_m)^2/(RFD_m - RFD_m2)
  #asymptotic value
  asy <- Diversity_profile(x,qs)
  #beta
  beta <- rep(0,length(qs))
  beta0plus <- which(asy != obs)
  beta[beta0plus] <- (Dn1[beta0plus]-obs[beta0plus])/(asy[beta0plus]-obs[beta0plus])
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
      RTD(ifi,n,m,qs) 
    }else if(m==n){
      obs
    }else{
      ETD(m,qs)
    }
  }
  sapply(m, Sub) %>% t() %>% as.vector()
}
TD.m.est_inc <- function(y, t_, qs){
  nT <- y[1]
  y <- y[-1]
  y <- y[y > 0]
  U <- sum(y)
  #xv_matrix = as.matrix(xv)
  iQi <- table(y);iQi <- cbind(i = as.numeric(names(iQi)),Qi = iQi)
  obs <- Diversity_profile_MLE.inc(c(nT,y),qs)
  RFD_m <- RTD_inc(iQi, nT, nT-1, qs)
  RFD_m2 <- RTD(iQi, nT, nT-2, qs)
  whincr <- which(RFD_m != RFD_m2)
  Dn1 <- obs; Dn1[whincr] <- obs + (obs - RFD_m)^2/(RFD_m - RFD_m2)
  asy <- Diversity_profile.inc(c(nT,y),qs)
  beta <- rep(0,length(qs))
  beta0plus <- which(asy != obs)
  beta[beta0plus] <- (Dn1[beta0plus]-obs[beta0plus])/(asy[beta0plus]-obs[beta0plus])
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
      RTD_inc(iQi,nT,m,qs) 
    }else if(m==nT){
      obs
    }else{
      ETD(m,qs)
    }
  }
  sapply(t_, Sub) %>% t() %>% as.vector()
}

AA.inc <- function(data){
  T = data[1]
  U <- sum(data[-1])
  data = data[-1]
  Yi = data[data!=0]
  Q1 <- sum(Yi==1)
  Q2 <- sum(Yi==2)
  if(Q2>0 & Q1>0){
    A <- 2*Q2/((T-1)*Q1+2*Q2)
  }
  else if(Q2==0 & Q1>1){
    A <- 2/((T-1)*(Q1-1)+2)
  }
  else{
    A <- 1
  }
  return(A)
}

