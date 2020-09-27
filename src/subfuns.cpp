#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double D1_2nd(double n, double f1, double A) {
  double q1 = 0;
  double h2 = 0;
  if(A==1||f1==0){
    h2 = 0;
  }else{
    for(int r = 1; r < n; r++){
      q1 = q1 + pow((1-A),r)/r;
    }
    h2 = (f1/n)*(pow(1-A,(-n+1)))*(-log(A)-q1);
  }
  return(h2);
}
// [[Rcpp::export]]
double Dq_2nd(double n, double f1, double A, double q) {
  double qq = 0;
  double ans = 0;
  if(A==1||f1==0){
    ans = 0;
  }else{
    for(int r = 0; r < n; r++){
      qq = qq + Rf_choose(q-1,r)*pow((A-1),r);
      //Rcpp::Rcout << "qq: " << qq << std::endl;
    }
    ans = (f1/n)*(pow(1-A,(-n+1)))*(pow(A,q-1)-qq);
  }
  return(ans);
}
// [[Rcpp::export]]
NumericVector Dq(NumericMatrix ifi, int n,NumericVector qs,double f1, double A){
  int nrows = ifi.nrow(), z = 0 , qlength = qs.length();
  double delta = 0.0;
  //NumericMatrix deltas(nrows,n);
  NumericMatrix ans_i(nrows,qlength);
  for(int i = 0; i<nrows; i++){
    z = ifi(i,0);
    for(int k = 0; k<=n-z; k++){
      delta = Rf_dhyper( 1, z, n - z, k+1, false )/(k+1);
      //deltas(i,k) = Rf_dhyper( 1, z, n - z, k+1, false )/(k+1);
      for(int i_q = 0; i_q<qlength; i_q++){
        ans_i(i,i_q) = ans_i(i,i_q) + ifi(i,1) * Rf_choose(k-qs[i_q],k)*delta;
      }
    }
  }
  NumericVector ans_1(qlength);
  for(int i_q = 0; i_q<qlength; i_q++){
    for(int i = 0; i<nrows; i++){
      ans_1(i_q) = ans_1(i_q) + ans_i(i,i_q);
    }
    ans_1(i_q) = ans_1(i_q) + Dq_2nd(n,f1,A,qs[i_q]);
    ans_1(i_q) = pow(ans_1(i_q),1/(1-qs[i_q]));
  }
  return(ans_1);
}

// [[Rcpp::export]]
NumericVector qDFUN(NumericVector q,NumericVector Xi,const int n){
  const int length = q.size();
  const int Sobs = Xi.size();
  NumericVector Q(length);
  NumericVector delta(n);
  NumericVector temp(Sobs);
  for(int k=0;k<=(n-1);k++){
    for(int i = 0;i<Sobs;i++){
      temp[i] = (Xi[i]/n)*exp(Rf_lchoose(n-Xi[i],k)-Rf_lchoose(n-1,k));
    }
    delta[k] = sum(temp);
  }
  
  for(int i=0;i<length;i++){
    float temp = 0;
    for(int k=0;k<=(n-1);k++){
      temp = temp + (Rf_choose(q[i]-1,k)*pow(-1,k)*delta[k]);
    }
    Q[i] = temp;
  }
  return Q;
}
// [[Rcpp::export]]
NumericVector qD_MLE(NumericVector q,NumericVector ai){
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
}

// [[Rcpp::export]]
NumericVector RTD(NumericMatrix x , int n  , double m , NumericVector q) {
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
}
// [[Rcpp::export]]
NumericVector RTD_inc(NumericMatrix y , int nT , double t_ , NumericVector q) {
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
}
