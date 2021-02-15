// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <vector>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]
double vectorSum(NumericVector x) {
  return std::accumulate(x.begin(), x.end(), 0.0, std::plus<double>());
}
// [[Rcpp::export]]
double MatrixSum(NumericMatrix A) {
  int n = A.nrow();
  NumericVector totvec(n);
  
  for(int i = 0; i < n; ++i ) {
    
    NumericVector vec(n);
    vec = A(i,_); 
    totvec[i] = vectorSum(vec);
  }
  double sumtot = 0;
  sumtot =   vectorSum(totvec);
  return sumtot;
}
 
  
  
  
//[[Rcpp::export]]
double makeTX2 (NumericMatrix X,
                  NumericMatrix tau,
                  double delta) {
    
    int n = X.ncol();
    int Q = tau.ncol();
    NumericMatrix QL(Q,Q);
    
    for(int q = 0; q < Q; ++q ) {
      for(int l = 0; l  <Q; ++l ) {
        
        if(q != l){
          NumericMatrix ij(n,n);
          for(int i = 0; i < n; ++i ) {
            for(int j = 0; j < n; ++j ) {
              
              if(i != j){
                double Meaned_Squared =  (X(i,j)-delta)*(X(i,j)-delta) ;
                ij(i,j) =  tau(i,q)*tau(j,l)* Meaned_Squared;
              }
              
            }
          }
          QL(q,l) =  MatrixSum(ij) ;
        }
      }  
    }
    return sum(QL);
    
  }


// [[Rcpp::export]]
NumericVector makeTX2_Q ( NumericMatrix tau ,
                  NumericMatrix X ,
                  double delta) {
  
  int n = tau.nrow();
  int Q = tau.ncol();
  NumericVector QQ(Q); 
  
  for(int q = 0; q < Q; ++q ) { 
    
    NumericMatrix ij(n,n);
    for(int i = 0; i < n; ++i ) {
      for(int j = 0; j < n; ++j ) {
        
        if(i != j){
          double Meaned_Squared =  (X(i,j)-delta)*(X(i,j)-delta) ;
          ij (i,j) =   tau(i,q)*tau(j,q)*Meaned_Squared   ;
        }
      }
    }
    QQ[q] = MatrixSum(ij) ;
  }
  return  QQ;
}


// [[Rcpp::export]]
double makeTX2_PQ ( NumericMatrix tau ,
                    NumericMatrix X ,
                    NumericVector P00,
                    double delta) {
  
  int n = tau.nrow();
  int Q = tau.ncol();
  NumericVector QQ(Q); 
  
  for(int q = 0; q < Q; ++q ) { 
    
    NumericMatrix ij(n,n);
    for(int i = 0; i < n; ++i ) {
      for(int j = 0; j < n; ++j ) {
        
        if(i != j){
          double Meaned_Squared =  (X(i,j)-delta)*(X(i,j)-delta) ;
          ij (i,j) =   tau(i,q)*tau(j,q)*Meaned_Squared* P00[q]  ;
        }
      }
    }
    QQ[q] = MatrixSum(ij) ;
  }
  return vectorSum(QQ);
}


//[[Rcpp::export]]
NumericVector MStep_SigmaX (NumericMatrix A,
                           NumericMatrix tau,
                           NumericVector Mu,
                           double xi ) {
  
  int n = A.ncol();
  int Q = tau.ncol();
  NumericVector QQ(Q);
  
  for(int q = 0; q < Q; ++q ) {
    NumericMatrix top(n,n);
    NumericMatrix bottom(n,n);
    double ratio=0;
    
    for(int i = 0; i < n; ++i ) {
      for(int j = 0; j < n; ++j ) {
        
        if(i != j ){
          
          double Signal = pow(A(i,j) - Mu[q], 2)  ;

          top(i,j) =  tau(i,q)*tau(j,q)* Signal ;
          bottom(i,j) =  tau(i,q)*tau(j,q)  ;
        }
        
      }
    }
    ratio =  MatrixSum(top)/MatrixSum(bottom) ;
    if (arma::is_finite(ratio)){
      QQ[q] =  ratio ;
    }else{
      QQ[q] =  0 ;
    }
  }
  return QQ;
}
 

 
 
 //[[Rcpp::export]]
 NumericVector MStep_SigmaXY    (NumericMatrix A1,
                                 NumericMatrix A2,
                                 NumericMatrix tau,
                                 NumericVector Mu1,
                                 NumericVector Mu2,
                                 double xi,
                                 double delta){
   
   int n = A1.ncol();
   int Q = tau.ncol();
   NumericVector QQ(Q);
   
   for(int q = 0; q < Q; ++q ) {
     NumericMatrix top(n,n);
     NumericMatrix bottom(n,n);
     double ratio=0;
     
     for(int i = 0; i < n; ++i ) {
       for(int j = 0; j < n; ++j ) {
         
         if(i != j){
           double Cov11 =  (A2(i,j) - Mu2[q]) * (A1(i,j) - Mu1[q])  ;
           top(i,j) =  tau(i,q)*tau(j,q)*Cov11 ;
           bottom(i,j) =  tau(i,q)*tau(j,q)    ;
         }
         
       }
     }
     
     ratio =  MatrixSum(top)/MatrixSum(bottom) ;
     if (arma::is_finite(ratio)){
       QQ[q] =  ratio ;
     }else{
       QQ[q] =  0 ;
     }
   }
   return QQ;
 }




// [[Rcpp::export]]
NumericVector makeT1_Q2 ( NumericMatrix tau ,
                          NumericVector P11) {
  
  int n = tau.nrow();
  int Q = tau.ncol();
  NumericVector QQ(Q); 
  
  for(int q = 0; q < Q; ++q ) { 
    
    NumericMatrix ij(n,n);
    for(int i = 0; i < n; ++i ) {
      for(int j = 0; j < n; ++j ) {
        
        if(i != j){
          ij (i,j) =   tau(i,q)*tau(j,q)*(  P11[q]) ;
        }
      }
    }
    QQ[q] = MatrixSum(ij) ;
  }
  return  QQ;
}

 