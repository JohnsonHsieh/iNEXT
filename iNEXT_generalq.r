TD.m.est_yhc = function(x, m, qs){
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

TD.m.est_inc_yhc <- function(y, t_, qs){
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

Diversity_profile <- function(x,q){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  r <- 1:(n-1)
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      A <- sum(x/n*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(-log(p1)-sum((1-p1)^r/r)))
      exp(A+B)
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      #ifelse(A==0,NA,A^(1/(1-q)))
      A^(1/(1-q))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      r <- 0:(n-1)
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      (A+B)^(1/(1-q))
    }
  }
  sapply(q, Sub)
}

Diversity_profile_MLE <- function(x,q){
  p <- x[x>0]/sum(x)
  Sub <- function(q){
    if(q==0) sum(p>0)
    else if(q==1) exp(-sum(p*log(p)))
    else exp(1/(1-q)*log(sum(p^q)))
  }
  sapply(q, Sub)
}

TD.m.est_yhc = function(x, m, qs){
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

cppFunction('NumericVector RTD(NumericMatrix x , int n  , double m , NumericVector q) {
  // x is x*i_vi; return length q
  int nrows = x.nrow();
  int qlength = q.length();
  NumericVector fhat(m);
  NumericVector out(qlength);
  //Rcout << "out in cpp: " << out << std::endl;
  for (int k = 0; k < m ; k++) {
    for (int i = 0; i < nrows; i++) {
      if (x(i,0) >= k+1 && x(i,0) <= n-m+k+1 )
      {
        fhat[k] += x(i,1)*exp(Rf_lchoose(x(i,0), k+1)+Rf_lchoose(n-x(i,0), m-k-1)-Rf_lchoose(n, m)) ;
      }
      else
      {
        fhat[k] += 0 ;
      }
    }
  }
  //Rcpp::Rcout << "hhat in cpp: " << hhat << std::endl;
  for (int j = 0; j < qlength; j++ ){  
    for(int k = 0; k < m; k++){
      if(q[j] == 0){
        out[j] = fhat[k] + out[j];
      }else if(q[j] == 1){
        //Rcout << "q1 in cpp: " <<log ( (k+1) )<< std::endl;
        out[j] = -( (k+1) / m ) * log ( (k+1) / m ) * fhat[k] + out[j];
      }else if(q[j] == 2){
        out[j] = pow( ( (k+1) / m ),2) * fhat[k] + out[j];
      }else{
        out[j] = pow( ( (k+1) / m ),q[j]) * fhat[k] + out[j];
      }
    }
  }
  //Rcout << "out in cpp: " << out << std::endl;
  for(int j = 0; j < qlength; j++ ){
    if(q[j] == 0){
      out[j] = out[j] ;
    }else if(q[j] == 1){
      out[j] = exp(out[j]);
    }else if(q[j] == 2){
      out[j] = 1 / out[j];
    }else{
      out[j] = pow( (out[j]) , 1/(1-q[j]) );
    }
  }
  return out;
}')
cppFunction('NumericVector RTD_inc(NumericMatrix y , int nT , double t_ , NumericVector q) {
  // x is x*i_vi; return length q
  int nrows = y.nrow();
  int qlength = q.length();
  NumericVector Qhat(t_);
  NumericVector out(qlength);
  //Rcout << "out in cpp: " << out << std::endl;
  for (int k = 0; k < t_ ; k++) {
    for (int i = 0; i < nrows; i++) {
      if (y(i,0) >= k+1 && y(i,0) <= nT-t_+k+1 )
      {
        Qhat[k] += y(i,1)*exp(Rf_lchoose(y(i,0), k+1)+Rf_lchoose(nT-y(i,0), t_-k-1)-Rf_lchoose(nT, t_)) ;
      }
      else
      {
        Qhat[k] += 0 ;
      }
    }
  }
  //Rcpp::Rcout << "hhat in cpp: " << hhat << std::endl;
  double U = 0;
  for (int i = 0; i < nrows; i++){
    U += y(i,1)*y(i,0);
  }
  double Ut_ = t_*U/nT;
  for (int j = 0; j < qlength; j++ ){  
    for(int k = 0; k < t_; k++){
      if(q[j] == 0){
        out[j] = Qhat[k] + out[j];
      }else if(q[j] == 1){
        //Rcout << "q1 in cpp: " <<log ( (k+1) )<< std::endl;
        out[j] = -( (k+1) / Ut_ ) * log ( (k+1) / Ut_ ) * Qhat[k] + out[j];
      }else if(q[j] == 2){
        out[j] = pow( ( (k+1) / Ut_ ),2) * Qhat[k] + out[j];
      }else{
        out[j] = pow( ( (k+1) / Ut_ ),q[j]) * Qhat[k] + out[j];
      }
    }
  }
  //Rcout << "out in cpp: " << out << std::endl;
  for(int j = 0; j < qlength; j++ ){
    if(q[j] == 0){
      out[j] = out[j] ;
    }else if(q[j] == 1){
      out[j] = exp(out[j]);
    }else if(q[j] == 2){
      out[j] = 1 / out[j];
    }else{
      out[j] = pow( (out[j]) , 1/(1-q[j]) );
    }
  }
  return out;
}')
