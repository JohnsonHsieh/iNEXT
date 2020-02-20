cppFunction(
  "NumericVector qD_MLE(NumericVector q,NumericVector ai){
    const int length = q.size();
    const int S = ai.size();
    NumericVector Q(length);
    NumericVector temp(S);
    for(int j = 0; j<length;j++){
    for(int i = 0 ; i<S;i++){
    temp[i] = pow(ai[i],q[j]);
    }
    Q[j] = pow(sum(temp),1/(1-q[j]));
    }
    return Q;
}")
Diversity_profile_MLE.inc <- function(data,q){
  Yi = data[-1]
  U = sum(Yi)
  Yi <- Yi[Yi!=0]
  ai <- Yi/U
  qD = qD_MLE(q,ai)
  qD[which(q==1)] <- exp(-sum(ai*log(ai)))
  return(qD)
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