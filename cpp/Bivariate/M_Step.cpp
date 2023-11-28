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
 
 
 
 
// [[Rcpp::export]]
double makeT0 ( NumericMatrix tau ) {
   
   int n = tau.nrow();
   int Q = tau.ncol();
   
   NumericMatrix QL(Q,Q); 
   
   for(int q = 0; q < Q; ++q ) {
     for(int l = 0; l < Q; ++l ) {
       
       if(q != l){
         NumericMatrix ij(n,n);
         for(int i = 0; i < n; ++i ) {
           for(int j = 0; j < n; ++j ) {
             
             if(i != j){
               ij (i,j) =   tau(i,q)*tau(j,l);
             }
           }
         }
         QL(q,l) = MatrixSum(ij) ;
       }
     }
   }
   return MatrixSum(QL);
 }

// [[Rcpp::export]]
double makeTX ( NumericMatrix tau ,  NumericMatrix X ) {
  
  int n = tau.nrow();
  int Q = tau.ncol();
  
  NumericMatrix QL(Q,Q); 
  
  for(int q = 0; q < Q; ++q ) {
    for(int l = 0; l < Q; ++l ) {
      
      NumericMatrix ij(n,n);
      
      if(q!=l){
        for(int i = 0; i < n; ++i ) {
          for(int j = 0; j < n; ++j ) {
     
             if(i != j){      
              ij (i,j) =   tau(i,q)*tau(j,l)*X(i,j);
             }
          }
        }
      }
      QL(q,l) = MatrixSum(ij) ;
    }
  }
  return MatrixSum(QL);
}


// [[Rcpp::export]]
double makeT0_Q ( NumericMatrix tau , NumericVector P00) {
  
  int n = tau.nrow();
  int Q = tau.ncol();
  NumericVector QQ(Q); 
  
  for(int q = 0; q < Q; ++q ) { 
    
    NumericMatrix ij(n,n);
    for(int i = 0; i < n; ++i ) {
      for(int j = 0; j < n; ++j ) {
        
        if(i != j){
          ij (i,j) =   tau(i,q)*tau(j,q)*(  P00[q]);
        }
      }
    }
    QQ[q] = MatrixSum(ij) ;
  }
  return vectorSum(QQ);
}

// [[Rcpp::export]]
NumericVector  makeT_Q ( NumericMatrix tau) {
  
  int n = tau.nrow();
  int Q = tau.ncol();
  NumericVector QQ(Q); 
  
  for(int q = 0; q < Q; ++q ) { 
    
    NumericMatrix ij(n,n);
    for(int i = 0; i < n; ++i ) {
      for(int j = 0; j < n; ++j ) {
        
        if(i != j){
          ij (i,j) =   tau(i,q)*tau(j,q) ; 
        }
      }
    }
    QQ[q] = MatrixSum(ij) ;
  }
  return QQ;
}



// [[Rcpp::export]]
NumericVector makeTX_Q ( NumericMatrix tau ,
                         NumericMatrix X ) {
  
  int n = tau.nrow();
  int Q = tau.ncol();
  NumericVector QQ(Q); 
  
  for(int q = 0; q < Q; ++q ) { 
    
    NumericMatrix ij(n,n);
    for(int i = 0; i < n; ++i ) {
      for(int j = 0; j < n; ++j ) {
        
        if( i != j){
          ij (i,j) =   tau(i,q)*tau(j,q)*X(i,j) ;
        }
      }
    }
    QQ[q] = MatrixSum(ij) ;
  }
  return  QQ;
}


// [[Rcpp::export]]
double makeTX_PQ ( NumericMatrix tau ,
                  NumericMatrix X ,
                  NumericVector P00) {
  
  int n = tau.nrow();
  int Q = tau.ncol();
  NumericVector QQ(Q); 
  
  for(int q = 0; q < Q; ++q ) { 
    
    NumericMatrix ij(n,n);
    for(int i = 0; i < n; ++i ) {
      for(int j = 0; j < n; ++j ) {
        
        if( i != j){
          ij (i,j) =   tau(i,q)*tau(j,q)*X(i,j)* ( P00[q]);
        }
      }
    }
    QQ[q] = MatrixSum(ij) ;
  }
  return vectorSum(QQ);
}



 // [[Rcpp::export]]
 NumericVector MStep_Mu_Q (NumericMatrix A,  
                         NumericMatrix tau ) {
   
   int n = A.ncol();
   int Q = tau.ncol();
   
   /* fill in w*/
   
   NumericVector QQ(Q);
   
   for(int q = 0; q < Q; ++q ) {
     for(int l = 0; l  <Q; ++l ) {
       
       NumericMatrix top(n,n);
       NumericMatrix bottom(n,n);
       double ratio=0;
       
       if(q==l){
         for(int i = 0; i < n; ++i ) {
           for(int j = 0; j < n; ++j ) {
             
             if(i != j){
               top(i,j) =  tau(i,q)*tau(j,l)* A(i,j) ;
               bottom(i,j) =  tau(i,q)*tau(j,l) ;
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
     }
   }
   
   return QQ;
 }

 