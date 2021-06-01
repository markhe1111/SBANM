# load packages
library(matlib)
library(DescTools)   
library(Brobdingnag)
library(igraph) 
library(MASS)
library(fpc)
library(rlist) 
library(randnet)

# load functions
Rcpp::sourceCpp('cpp/Trivariate/E_Step.cpp')
Rcpp::sourceCpp('cpp/Trivariate/M_Step.cpp')
Rcpp::sourceCpp('cpp/Trivariate/M_Step_b.cpp')

source('SBANM/Utilities/general_Utils.R')
source('SBANM/Utilities/Triv_Utils.R')
source('SBANM/Trivariate.R')

#INPUT: 
# A1,A2,A3 are weighted adjacency matrices
# Q is selection of number of blocks 
 # TOL is tolerance for cutting off 


trivSBANM = function(Q, A1,  A2, A3, TOL =1e-12 ){
  
  
  A_sum = Reduce('+' , list(A1, A2, A3)) 
  sphericalspectra_tot = reg.SSP(A_sum,Q) 
  Psi = (Q-1)/Q
  
  n_test = ncol(A1)
  for(initializing_W_Test in 1:1){
    ELBOs = list()
    bigbigtau_list = list()
    bigbigP._list = list()
    n = nrow(A1) 
    Alpha_list = list()
    bigtau_list = list();  tau_list = list() ;   alpha_list = list();
    taus = list();  diffs = list() ;   P_list = list()
    mu_list = list();  D.log = list()
    
    initial_1=sapply(sphericalspectra_tot$cluster, function(x)  ones (x,Q) )
   
    tau.01 =  t(initial_1) ;  
    tau.02 =   t(  sapply(1:n_test,  function(x) rep(1/Q,Q)) )
    tau.0 = tau.01 + tau.02  ;  
    tau. = tau.0 /  rowSums(tau.0)
    
    tau =tau.
    
    P. = cbind(rep( 1- 1/Q,Q), rep( 1/Q ,Q))
    
    P.11 = P.[,1];  P.00 = P.[,2]
    alpha = colMeans(tau)
    
    TX = makeTX(tau, A1)
    TY = makeTX(tau, A2)
    TZ = makeTX(tau, A3)
    
    T0 = makeT0(tau )
    TX_Q = makeTX_PQ (tau, A1, P.00 )
    TY_Q = makeTX_PQ (tau, A2, P.00)
    TZ_Q = makeTX_PQ (tau, A3, P.00)
    
    T0_Q = makeT0_Q(tau,  P.00) 
    
    xi1 =   (Psi)* nan_omit(TX / T0)    + (1-Psi)* nan_omit( TX_Q /  T0_Q )
    xi2 =   (Psi)* nan_omit(TY / T0)    + (1-Psi)* nan_omit( TY_Q /  T0_Q )
    xi3 =   Psi* nan_omit(TZ / T0)    + (1-Psi)* nan_omit( TZ_Q /  T0_Q )
    
    Xis = c(xi1,xi2,xi3)
    
    Mu1 = MStep_Mu_Q(A1,tau) * P.11 + (1-P.11)* xi1
    Mu2 = MStep_Mu_Q(A2,tau  )* P.11+ (1-P.11)* xi1
    Mu3 = MStep_Mu_Q(A3,tau  )* P.11+ (1-P.11)* xi1
    
    
    TX2 =  makeTX2(A1,tau,xi1)
    TY2 =  makeTX2(A2,tau,xi2)
    TZ2 =  makeTX2(A3,tau,xi3)
    
    # PQ is in alt
    
    TX2Q = makeTX2_PQ  (tau,A1, P00 = P.00   ,xi1)
    TY2Q = makeTX2_PQ  (tau,A2, P00 = P.00   ,xi2)
    TZ2Q = makeTX2_PQ  (tau,A3, P00 = P.00   ,xi3)
    
    Sigmaxi1  =    (Psi)* nan_omit(TX2 /T0)  +  (1-Psi)* nan_omit(TX2Q /  T0_Q )
    Sigmaxi2  =   (Psi)*  nan_omit(TY2 /T0)  +  (1-Psi)*nan_omit(TY2Q /  T0_Q )
    Sigmaxi3  =   (Psi)*  nan_omit(TZ2 /T0 ) +  (1-Psi)*nan_omit(TZ2Q /  T0_Q )
    
    sigmaXis = c(sqrt(Sigmaxi1),sqrt(Sigmaxi2),sqrt(Sigmaxi3))
    
    SigmaX = MStep_SigmaX ( A1,tau,Mu1, xi1 ) * P.11 + Sigmaxi1*P.00
    SigmaY = MStep_SigmaX ( A2,tau,Mu2, xi2)* P.11 + Sigmaxi2*P.00
    SigmaZ = MStep_SigmaX ( A3,tau,Mu3, xi3)* P.11 + Sigmaxi3*P.00
    
    SigmaXY = MStep_SigmaXY  (A1,A2,tau,Mu1,Mu2 ,xi1,xi2 )*P.11
    SigmaYZ = MStep_SigmaXY  (A2,A3,tau,Mu2,Mu3 ,xi2,xi3 )*P.11
    SigmaXZ = MStep_SigmaXY (A1,A3,tau,Mu1,Mu3 ,xi1,xi3 )*P.11
    
    rhos12 = SigmaXY  / sqrt(SigmaX * SigmaY)
    rhos23 = SigmaYZ  / sqrt(SigmaY * SigmaZ)
    rhos13 = SigmaXZ  / sqrt(SigmaX * SigmaZ)
    rhos12[is.na(rhos12)] = 0
    rhos23[is.na(rhos23)] = 0
    rhos13[is.na(rhos13)] = 0
    
  }
  options(warn=2)
  
  for(h in  1 :200 ){
    
    tau_prev = tau
    delta_t =  1- 1/(h+1) 
    StochSize =   min(n, round(  100 + n *  delta_t^2)  )
    m= sort(sample(n,  StochSize))
    tau.m = tau[m,]
    
    QQs = sapply( 1:Q,   function(q) get_EN_each_q(q, A1 ,   A2,  A3 , tau ,  
                                                   Mu1,  Mu2,  Mu3,  Xis, 
                                                   SigmaX, SigmaY, SigmaZ,
                                                   SigmaXY, SigmaYZ,  SigmaXZ, sigmaXis))
    
    
    
    
    Nq_inv1 =  as.brob( exp( QQs + log( (1-Psi)/  Psi)) )
    Nq_inv2 =   ( 1 + 1/ (Nq_inv1) )
    
    N_n = 1/ Nq_inv2
    norm_N =  nan_omit ( as.numeric( N_n/ sum(N_n) ))
    
    if(sum(norm_N)==0){
      norm_N = rep(1/Q,Q)
    }
    
    P.. =   cbind(1- norm_N ,  norm_N)
    P. = P.. * delta_t + P. * (1-delta_t)
    P.11 = P.[,1] ; P.00 = P.[,2]
    
    taus._m = list()
    for(t in 1:30){
      t12 = E_tau_Triv_Mat (  A1 [m,m],   A2[m,m],  A3[m,m], tau[m,] ,  
                              Mu1,  Mu2,  Mu3, SigmaX, SigmaY, SigmaZ,
                              SigmaXY, SigmaYZ,  SigmaXZ,
                              Xis ,   sigmaXis,     P1 =P.11)
      t12[is.nan(t12)] = -Inf
      
      tsum1 = sapply(1:Q, function(q) t12[,q] + log(alpha)[q]  ) 
      tau.m = NormalizeTauBrob(tsum1 )
      tau.m[is.nan(tau.m)] = 1/Q 
      taus._m[[t]] = tau.m
      
      if(t>1) {
        diffs[[t]] = sum(abs( tau.m -taus._m[[t-1]]))
        
        if(   diffs[[t]]  > 100 ) break 
        if(   diffs[[t]]  < 1e-6  ) break
        if(t>3){
          if(  abs( (diffs[[t]] - diffs[[t-1]]  )/ diffs[[t-1]] )   < 1e-6 ) break
        }
        print (paste('tau diff - micro iter', diffs[[t]]))
      }
    }
    
    tau[m,]  = tau_prev[m,] * (1-delta_t) + tau.m * delta_t
    alpha = colMeans(tau)
    
    TX = makeTX(tau, A1)
    TY = makeTX(tau, A2)
    TZ = makeTX(tau, A3)
    
    T0 = makeT0(tau )
    TX_Q = makeTX_PQ (tau, A1, P.00 )
    TY_Q = makeTX_PQ (tau, A2, P.00)
    TZ_Q = makeTX_PQ (tau, A3, P.00)
    
    T0_Q = makeT0_Q(tau,  P.00) 
    
    xi1 =   (Psi)* nan_omit(TX / T0)    + (1-Psi)* nan_omit( TX_Q /  T0_Q )
    xi2 =   (Psi)* nan_omit(TY / T0)    + (1-Psi)* nan_omit( TY_Q /  T0_Q )
    xi3 =     Psi* nan_omit(TZ / T0)    + (1-Psi)* nan_omit( TZ_Q /  T0_Q )
    
    
    Xis = c(xi1,xi2,xi3)
    
    Mu1 = MStep_Mu_Q(A1,tau) * P.11 +(1-P.11)* xi1
    Mu2 = MStep_Mu_Q(A2,tau )* P.11 + (1-P.11)* xi2
    Mu3 = MStep_Mu_Q(A3,tau )* P.11 + (1-P.11)* xi3
    
    Mu1[P.11==0] = xi1
    Mu2[P.11==0] = xi2
    Mu3[P.11==0] = xi3
    
    TX2 =  makeTX2(A1,tau,xi1)
    TY2 =  makeTX2(A2,tau,xi2)
    TZ2 =  makeTX2(A3,tau,xi3)
    
    TX2Q = makeTX2_PQ  (tau,A1,P.00, xi1)
    TY2Q = makeTX2_PQ  (tau,A2,P.00, xi2)
    TZ2Q = makeTX2_PQ  (tau,A3,P.00, xi3) 
    
    Sigmaxi1  =   (Psi)* nan_omit(TX2 /T0)  +   (1-Psi)* nan_omit(  (TX2Q) /  T0_Q )
    Sigmaxi2  =   (Psi)*  nan_omit(TY2 /T0)  +  (1-Psi)*nan_omit(   (TY2Q) /  T0_Q )
    Sigmaxi3  =   (Psi)*  nan_omit(TZ2 /T0 ) +  (1-Psi)*nan_omit(   (TZ2Q) /  T0_Q )
    
    sigmaXis = c(sqrt(Sigmaxi1),sqrt(Sigmaxi2),sqrt(Sigmaxi3))
    
    SigmaX = MStep_SigmaX ( A1,tau,Mu1, xi1 )*P.11 + Sigmaxi1*P.00
    SigmaY = MStep_SigmaX ( A2,tau,Mu2, xi2)* P.11 + Sigmaxi2*P.00
    SigmaZ = MStep_SigmaX ( A3,tau,Mu3, xi3)* P.11 + Sigmaxi3*P.00
    
    SigmaXY = MStep_SigmaXY  (A1,A2,tau,Mu1,Mu2 ,xi1,xi2 )*P.11
    SigmaYZ = MStep_SigmaXY  (A2,A3,tau,Mu2,Mu3 ,xi2,xi3 )*P.11
    SigmaXZ = MStep_SigmaXY (A1,A3,tau,Mu1,Mu3 ,xi1,xi3 )*P.11
    
    SigmaX[P.11==0] = Sigmaxi1
    SigmaY[P.11==0] = Sigmaxi2
    SigmaZ[P.11==0] = Sigmaxi3
    
    rhos12 =  ( nan_omit( SigmaXY  / sqrt(SigmaX * SigmaY)))
    rhos23 =  ( nan_omit(SigmaYZ  / sqrt(SigmaY * SigmaZ)) )
    rhos13 =   ( nan_omit(SigmaXZ  / sqrt(SigmaX * SigmaZ)) )
    
    
    rhos_temp = cbind(rhos12, rhos23,rhos13 )
    rhos12 = rowMax(rhos_temp) 
    rhos13 = rowMax(rhos_temp) 
    rhos23 = rowMax(rhos_temp) 
    rhosmax = cbind(rhos12, rhos23,rhos13 )
    rhosmax[rhosmax<0] = 0
    
    CrossList = make_variances_Trivariate  (Q , SigmaX, SigmaY, SigmaZ,rhosmax)
    
    SigmaXY = sapply( CrossList, function(x)   x[1,2])
    SigmaYZ = sapply( CrossList, function(x)   x[2,3])
    SigmaXZ = sapply( CrossList, function(x)   x[1,3])
    
    Alpha_list[[h]] = alpha
    bigtau_list[[h]] = tau
    P_list[[h]] = P. 
    
    rhos = cbind(rhos12, rhos23,rhos13 )
    
    if(h>1){
      print (list('all tau', sum(abs( bigtau_list[[h]] -  bigtau_list[[h-1]]   )) , 
                  tautable(tau),P.,round(rhos,2) ))
      if (sum(abs( bigtau_list[[h]] -  bigtau_list[[h-1]] )) < 1e-6 ) break
    }else{
      print(paste('end first round', tautable(tau)))
    }
  }  
  
  RunList = list(    tau=tau, Mu1=Mu1, Mu2=Mu2,Mu3=Mu3, P.11=P.11, rhos=rhos ,
                     SigmaXY=SigmaXY, SigmaYZ=SigmaYZ, SigmaXZ=SigmaXZ,
                     SigmaX =SigmaX, SigmaY=SigmaY, SigmaZ=SigmaZ,
                     Xis = Xis, sigmaXis=sigmaXis,
                     Alpha_list=Alpha_list, bigtau_list=bigtau_list ,P_list=P_list )
  return(RunList)
}



