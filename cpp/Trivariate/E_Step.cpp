// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins("cpp11")]]
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <vector>
#include <iostream>
// #include <bigmemory/MatrixAccessor.hpp>


using namespace Rcpp;

static const double pi = 3.1415926535897932384626; 

// https://github.com/petewerner/misc/wiki/RcppArmadillo-cheatsheet


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
 
 // [[ Rcpp :: export ()]]
 double Mprod (arma::mat x,
                arma::vec y) {
   arma::mat ip = x.t() * y ;
   return(ip(0)) ;
 }
 

 // [[Rcpp::export]]
NumericVector VecDiff(NumericVector ys, NumericVector xs) {
   int n = ys.size();
   NumericVector out(n);
   
   for(int i = 0; i < n; ++i) {
     out[i] =  ys[i] - xs[i];
   }
   return out;
 }
 
 
 
 
// [[Rcpp::export]]
double getTrivSignal( NumericVector X,  NumericVector Mu, 
                      arma::mat   Sigma){
  // arma::mat  Sigma
  arma::mat InvSigma = arma::pinv(Sigma) ;
  arma::vec Xmu = VecDiff(X,Mu) ;
  
  arma::mat InvSigX1 =  InvSigma * Xmu  ;
  double InvSigX2 = Mprod(InvSigX1 , Xmu ); 
  
  double logterm =  pow(2*pi,3) * arma::det(Sigma)    ;
  double freqtot = - .5* InvSigX2 - .5*log(logterm);
   
  return freqtot;
} 





// [[Rcpp::export]]
arma::mat makeCovSig(double q, 
                            NumericVector SigmaX, NumericVector SigmaY , NumericVector SigmaZ,
                            NumericVector  SigmaXY, NumericVector SigmaYZ ,  NumericVector SigmaXZ){
   
  
  arma::mat CovSig (3,3);
  
  CovSig(0,0) = SigmaX[q];
  CovSig(1,1) = SigmaY[q];
  CovSig(2,2) = SigmaZ[q];
  
  CovSig(1,0) =   SigmaXY[q];
  CovSig(0,1) =   SigmaXY[q];
  
  CovSig(1,2) =   SigmaYZ[q];
  CovSig(2,1) =   SigmaYZ[q];
  
  CovSig(2,0) =   SigmaXZ[q];
  CovSig(0,2) =   SigmaXZ[q];

  return CovSig;
  
}

 


// [[Rcpp::export]]
NumericVector E_tau_Triv_Mat    (NumericMatrix A1, 
                                 NumericMatrix A2,  
                                 NumericMatrix A3,  
                                 NumericMatrix tau,  
                                 NumericVector Mu1,
                                 NumericVector Mu2,
                                 NumericVector Mu3,
                                 NumericVector SigmaX,
                                 NumericVector SigmaY,
                                 NumericVector SigmaZ,
                                 NumericVector SigmaXY,
                                 NumericVector SigmaYZ,
                                 NumericVector SigmaXZ,
                                 NumericVector Xis,
                                 NumericVector sigmaXis, 
                                 NumericVector P1 ){
  
  int n = tau.nrow();
  int Q = tau.ncol();
  
  /* fill in w*/
  NumericMatrix JL (n,Q);
  
  for(int i = 0; i < n; ++i ) {
    for(int q = 0; q < Q; ++q ) {
      
      NumericMatrix IQ (n,Q);  
      double IQ_sum;
      
      for(int l = 0; l < Q; ++l ) {
        for(int j = 0; j < n; ++j ) {
          
          if(i != j){ 
            
            double IQjl = 0; 
            NumericVector X = NumericVector::create(  A1(i,j), A2(i,j), A3(i,j) );
            
            if(q==l){   
              
              NumericVector Mus = NumericVector::create(  Mu1[q], Mu2[q], Mu3[q]);
              arma::mat CovMat =   makeCovSig( q, SigmaX,SigmaY , SigmaZ,
                                                        SigmaXY, SigmaYZ ,  SigmaXZ);

              arma::mat  CovNoi (3,3);
              CovNoi(0,0) = pow(sigmaXis[0],2) ;
              CovNoi(1,1) = pow(sigmaXis[1],2) ;
              CovNoi(2,2) = pow(sigmaXis[2],2) ;  
              
              double SignalTerm = getTrivSignal(X,Mus, CovMat)  ;
              double NoiseTerm = getTrivSignal(X,Xis, CovNoi)  ;
                
              // log term:
              IQjl =  SignalTerm *P1[q] + NoiseTerm*(1-P1[q])  ; 
              
            }else{ 
              
              arma::mat CovNoi (3,3);
              CovNoi(0,0) = pow(sigmaXis[0],2) ;
              CovNoi(1,1) = pow(sigmaXis[1],2) ;
              CovNoi(2,2) = pow(sigmaXis[2],2) ;  
              
              double NoiseTerm = getTrivSignal(X,Xis, CovNoi)  ;
              IQjl =    NoiseTerm ;
            } 
            
            IQ(j,l) =   tau(j,l) * IQjl;
          }
        }
      }
      IQ_sum = MatrixSum(IQ)  ;
      JL(i,q) =   IQ_sum;
      
    }
  }
  
  return JL;
}




// [[Rcpp::export]]
NumericVector E_P   (NumericMatrix A1, 
                     NumericMatrix A2,  
                     NumericMatrix A3,  
                     NumericMatrix tau,  
                     NumericVector Mu1,
                     NumericVector Mu2,
                     NumericVector Mu3,
                     NumericVector SigmaX,
                     NumericVector SigmaY,
                     NumericVector SigmaZ,
                     NumericVector SigmaXY,
                     NumericVector SigmaYZ,
                     NumericVector SigmaXZ,
                     NumericVector Xis,
                     NumericVector sigmaXis ){
  int n = tau.nrow();
  int Q = tau.ncol();
  
  NumericVector QQ(Q); 
  
  for(int q = 0; q < Q; ++q ) { 
    
    NumericMatrix Pq_ij (n,n);
    
    NumericVector Mus = NumericVector::create(  Mu1[q], Mu2[q], Mu3[q]);
    arma::mat CovMat =   makeCovSig( q, SigmaX,SigmaY , SigmaZ, SigmaXY, SigmaYZ ,  SigmaXZ);
    arma::mat  CovNoi (3,3);
    CovNoi(0,0) = pow(sigmaXis[0],2) ;
    CovNoi(1,1) = pow(sigmaXis[1],2) ;
    CovNoi(2,2) = pow(sigmaXis[2],2) ;  
    
    for(int i = 0; i < n; ++i ) {
      for(int j = 0; j < n; ++j ) {
        
        if( i != j){
          NumericVector X = NumericVector::create(  A1(i,j), A2(i,j), A3(i,j) );
          double SignalTerm = getTrivSignal(X,Mus, CovMat)  ;
          double NoiseTerm  = getTrivSignal(X,Xis, CovNoi)  ;
          Pq_ij (i,j) =  tau(i,q)*tau(j,q)/2* (SignalTerm  -  NoiseTerm )    ; 
        }
      }
    }
    QQ[q] = MatrixSum(Pq_ij)   ;
  }
  return QQ;
}

  
  
  // [[Rcpp::export]]
  NumericVector E_Nq   (NumericMatrix A1, 
                       NumericMatrix A2,  
                       NumericMatrix A3,  
                       NumericMatrix tau,  
                       NumericVector Mu1,
                       NumericVector Mu2,
                       NumericVector Mu3,
                       NumericVector SigmaX,
                       NumericVector SigmaY,
                       NumericVector SigmaZ,
                       NumericVector SigmaXY,
                       NumericVector SigmaYZ,
                       NumericVector SigmaXZ,
                       NumericVector Xis,
                       NumericVector sigmaXis ){
    int n = tau.nrow();
    int Q = tau.ncol();
    
    NumericVector QQ(Q); 
    
    for(int q = 0; q < Q; ++q ) { 
      
      NumericMatrix Pq_ij (n,n);
      
      NumericVector Mus = NumericVector::create(  Mu1[q], Mu2[q], Mu3[q]);
      arma::mat CovMat =   makeCovSig( q, SigmaX,SigmaY , SigmaZ, SigmaXY, SigmaYZ ,  SigmaXZ);
      arma::mat  CovNoi (3,3);
      CovNoi(0,0) = pow(sigmaXis[0],2) ;
      CovNoi(1,1) = pow(sigmaXis[1],2) ;
      CovNoi(2,2) = pow(sigmaXis[2],2) ;  
      
      for(int i = 0; i < n; ++i ) {
        for(int j = 0; j < n; ++j ) {
          
          if( i != j){
            NumericVector X = NumericVector::create(  A1(i,j), A2(i,j), A3(i,j) );
            double SignalTerm = getTrivSignal(X,Mus, CovMat)  ;
            double NoiseTerm  = getTrivSignal(X,Xis, CovNoi)  ;
            Pq_ij (i,j) =  tau(i,q)*tau(j,q)/2* (NoiseTerm -SignalTerm  )    ; 
          }
        }
      }
      QQ[q] = MatrixSum(Pq_ij)   ;
    }
    return QQ;
  }




// [[Rcpp::export]]
double getTrivSignal_Inv( NumericVector X,  
                          NumericVector Mu, 
                          arma::mat   Sigma,
                          arma::mat   InvSigma){ 
  arma::vec Xmu = VecDiff(X,Mu) ;
  arma::mat InvSigX1 =  InvSigma * Xmu  ;
  double InvSigX2 = Mprod(InvSigX1 , Xmu ); 
  
  double logterm =  pow(2*pi,3) * arma::det(Sigma)    ;
  double freqtot = - .5* InvSigX2 - .5*log(logterm);
  return freqtot;
} 


// [[Rcpp::export]]
double E_Nq_Simp_q   (NumericMatrix A1, 
                      NumericMatrix A2,  
                      NumericMatrix A3,  
                      NumericVector tau_q,  
                      NumericVector Mus_q,
                      arma::mat CovMat_q,
                      arma::mat  InvCovMat_q,
                      NumericVector Xis,
                      arma::mat CovNoi,
                      arma::mat  InvCovNoi){
  int n = A2.nrow();

  double Qq; 
  NumericMatrix Pq_ij (n,n);
  
  for(int i = 0; i < n; ++i ) {
    for(int j = 0; j < n; ++j ) {
      
      if( i != j){
        NumericVector X = NumericVector::create(  A1(i,j), A2(i,j), A3(i,j) );
        double SignalTerm = getTrivSignal_Inv(X,Mus_q, CovMat_q, InvCovMat_q)  ;
        double NoiseTerm  = getTrivSignal_Inv(X,Xis,   CovNoi  , InvCovNoi)  ;
        Pq_ij (i,j) =  tau_q[i]*tau_q(j)/2* (NoiseTerm -SignalTerm  )    ; 
      }
    }
  }
  Qq = MatrixSum(Pq_ij)   ;
  return Qq;
}



// [[Rcpp::export]]
double E_Pq_Simp_q   (NumericMatrix A1, 
                      NumericMatrix A2,  
                      NumericMatrix A3,  
                      NumericVector tau_q,  
                      NumericVector Mus_q,
                      arma::mat CovMat_q,
                      arma::mat  InvCovMat_q,
                      NumericVector Xis,
                      arma::mat CovNoi,
                      arma::mat  InvCovNoi){
  int n = A2.nrow();
  
  double Qq; 
  NumericMatrix Pq_ij (n,n);
  
  for(int i = 0; i < n; ++i ) {
    for(int j = 0; j < n; ++j ) {
      
      if( i != j){
        NumericVector X = NumericVector::create(  A1(i,j), A2(i,j), A3(i,j) );
        double SignalTerm = getTrivSignal_Inv(X,Mus_q, CovMat_q, InvCovMat_q)  ;
        double NoiseTerm  = getTrivSignal_Inv(X,Xis,   CovNoi  , InvCovNoi)  ;
        Pq_ij (i,j) =  tau_q[i]*tau_q(j)/2* (SignalTerm-NoiseTerm   )    ; 
      }
    }
  }
  Qq = MatrixSum(Pq_ij)   ;
  return Qq;
}







