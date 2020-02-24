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
  nT = data[1]
  Yi = data[-1]
  Yi <- Yi[Yi!=0]
  U <- sum(Yi)
  Q1 <- sum(Yi==1)
  Q2 <- sum(Yi==2)
  Sobs <- length(Yi)
  A <- AA.inc(data)
  Q0hat <- ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2)
  B <- sapply(q,function(q) ifelse(A==1,0,(Q1/nT)*(1-A)^(-nT+1)*(A^(q-1)-sum(sapply(c(0:(nT-1)),function(r) choose(q-1,r)*(A-1)^r)))))
  qD <- (U/nT)^(q/(q-1))*(qDFUN(q,Yi,nT) + B)^(1/(1-q))
  qD[which(q==0)] = Sobs+Q0hat
  yi <- Yi[Yi>=1 & Yi<=(nT-1)]
  delta <- function(i){
    (yi[i]/nT)*sum(1/c(yi[i]:(nT-1)))
  }
  if(sum(q %in% 1)>0){
    C_ <- ifelse(A==1,0,(Q1/nT)*(1-A)^(-nT+1)*(-log(A)-sum(sapply(c(1:(nT-1)),function(r) (1-A)^r/r))))
    qD[which(q==1)] <- exp((nT/U)*( sum(sapply(c(1:length(yi)),function(i) delta(i))) + C_)+log(U/nT))
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

AA.inc <- function(data){
  nT = data[1]
  U <- sum(data[-1])
  data = data[-1]
  Yi = data[data!=0]
  Q1 <- sum(Yi==1)
  Q2 <- sum(Yi==2)
  if(Q2>0 & Q1>0){
    A <- 2*Q2/((nT-1)*Q1+2*Q2)
  }
  else if(Q2==0 & Q1>1){
    A <- 2/((nT-1)*(Q1-1)+2)
  }
  else{
    A <- 1
  }
  return(A)
}

