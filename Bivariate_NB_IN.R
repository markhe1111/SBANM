
library(matlib)
library(DescTools)   
library(Brobdingnag)
library(igraph) 
library(MASS)
library(fpc)
library(rlist) 
library(randnet)

Rcpp::sourceCpp('cpp/Bivariate/E_Biv_IN_NB.cpp')
Rcpp::sourceCpp('cpp/Bivariate/MStepBiv.cpp')
Rcpp::sourceCpp('cpp/Bivariate/biv_ELBO.cpp')

source('Utilities/general_Utils.R')
source('Utilities/Biv_Utils.R')

#INPUT: 
# A1,A2  are weighted adjacency matrices
# Q is selection of number of blocks 
# init is choice of initialization: 'spec' for Spectral clustering, and 'random' for random uniform assignment of memberships
biv_EN_each_q <-function(q, A1 ,  A2 , tau ,  
                         Mu1,  Mu2, xiBs,
                         SigmaX, SigmaY,
                         SigmaXY,  sigmaXiBs){
  
  CovMat_q = makeCovSig( q-1, SigmaX,SigmaY  , SigmaXY )  
  covNoi = diag( sigmaXiBs^2) 
  Mus_q = cbind(Mu1,Mu2)[q,]
  Freq_q =   E_Nq_Simp_Biv  ( A1,A2,   tau[,(q)],
                              Mus_q, CovMat_q,  ginv(CovMat_q),  
                              xiBs,   covNoi,    ginv(covNoi))
  Freq_q
}



bivSBANM = function(Q, A1,  A2, init="mix",
                    max_iter = 200, StochSize.a = 200 , TOL =1e-5 ){
  
  n = ncol(A1)
  A_sum = Reduce('+' , list(A1, A2 ))
  sphericalspectra_tot = reg.SSP(A_sum,Q)
  Psi = (Q-1)/Q
  
  for(initializing_W_Test in 1:1){
    ELBOs = list()
    bigbigtau_list = list()
    bigbigP._list = list()
    n = nrow(A1)
    Alpha_list = list()
    bigtau_list = list();  tau_list = list() ;   alpha_list = list();
    taus = list();  diffs = list() ;   P_list = list()
    mu_list = list();  D.log = list()
    
    if(init=='random'){
      initial = sapply(1:Q, function(x) runif(n, min=0, max=10) )
      tau = initial/rowSums(initial)
      
    } else if (init=='mix') {
      
      initial_1=sapply(sphericalspectra_tot$cluster, function(x)  ones (x,Q) )
      #initial_unif = sapply(1:Q, function(x) runif(n, min=0, max=10) )
      initial_2 = rep( 1/Q,Q ) ; #initial_2 = t(initial_unif/rowSums(initial_unif))
      initial = (initial_1 + initial_2) / colSums(initial_1+initial_2)
      tau = t(initial)
      
    }else if (init=='mix2') {
      
      initial_1=sapply(sphericalspectra_tot$cluster, function(x)  ones (x,Q) )
      initial_unif = sapply(1:Q, function(x) runif(n, min=0, max=10) )
      
      initial_2 = rep( 1/Q,Q ) ; 
      initial_3 = t(initial_unif/rowSums(initial_unif))
      
      
      initial = (initial_1 + initial_2 + initial_3) / colSums(initial_1+initial_2+initial_3)
      tau = t(initial)
    }
    else if (init=='plugin'){
      initial = init_plug / rowSums(init_plug)  
      
    }else{
      initial= sapply(sphericalspectra_tot$cluster, function(x)  ones (x,Q) )
      tau = t(initial  )
    }
    P. = cbind(rep( (Q-1)/Q  ,Q), rep(1/Q,Q))
    P.11 = P.[,1];  P.00 = P.[,2]
    alpha = colMeans(tau)
    
    TX = makeTX(tau, A1);  TY = makeTX(tau, A2)
    T0 = makeT0(tau )
    TX_Q = makeTX_PQ (tau, A1, P.00 )
    TY_Q = makeTX_PQ (tau, A2, P.00)
    
    T0_Q = makeT0_Q(tau,  P.00)
    
    xi1 =    nan_omit(TX / T0)   ;xi2 =    nan_omit(TY / T0)    ; xis = c(xi1,xi2)
    xiB1 =     nan_omit( TX_Q /  T0_Q ); xiB2 =     nan_omit( TY_Q /  T0_Q );  xiBs = c(xiB1,xiB2)
    Mu1 = MStep_Mu_Q(A1,tau) * P.11 + (1-P.11)* xi1
    Mu2 = MStep_Mu_Q(A2,tau  )* P.11+ (1-P.11)* xi2
    
    TX2 =  makeTX2(A1,tau,xi1);TY2 =  makeTX2(A2,tau,xi2)
    TX2Q = makeTX2_PQ  (tau,A1, P00 = P.00   ,xi1);  TY2Q = makeTX2_PQ  (tau,A2, P00 = P.00   ,xi2)
    
    Sigmaxi1  =   nan_omit(TX2 /T0) ; Sigmaxi2  =   nan_omit(TY2 /T0)   
    sigmaXis = c(sqrt(Sigmaxi1),sqrt(Sigmaxi2) )
    
    SigmaxiB1  =     nan_omit(TX2Q /  T0_Q ) ; SigmaxiB2  =     nan_omit(TY2Q /  T0_Q )
    sigmaXiBs = c(sqrt(SigmaxiB1),sqrt(SigmaxiB2) )
    
    SigmaX = MStep_SigmaX ( A1,tau,Mu1, xi1 ) * P.11 + Sigmaxi1*P.00
    SigmaY = MStep_SigmaX ( A2,tau,Mu2, xi2)* P.11 + Sigmaxi2*P.00
    
    SigmaXY = MStep_SigmaXY  (A1,A2,tau,Mu1,Mu2 ,xi1,xi2 )*P.11
    
    rhos12 = SigmaXY  / sqrt(SigmaX * SigmaY)
    rhos12[is.na(rhos12)] = 0
  }
  
  
  tt = proc.time()
  for(h in  1  :max_iter ){
    
    tau_prev = tau
    delta_t =  1- 1/(h+1)
    StochSize =   min(n, round(  500 + n *  delta_t^2)  )
    m= sort(sample(n,  StochSize))
    tau.m = tau[m,]
    
    QQs = sapply( 1:Q,   function(q) biv_EN_each_q(q, A1 ,  A2 , tau ,  
                                                   Mu1,  Mu2, xiBs,
                                                   SigmaX, SigmaY,
                                                   SigmaXY,  sigmaXiBs))
    Nq_inv1 =  as.brob( exp( QQs + log( (1-Psi)/  Psi)) )
    #Nq_inv1 =  as.brob( exp( QQs  ))
    
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
      t12 = E_tau_Biv_Mat (  A1 [m,m],   A2[m,m],   tau[m,] ,  
                             Mu1,  Mu2, SigmaX, SigmaY,
                             SigmaXY,
                             xis ,   sigmaXis,  
                             xiBs, sigmaXiBs,
                             P1 =P.11)
      t12[is.nan(t12)] = -Inf
      
      tsum1 = sapply(1:Q, function(q) t12[,q] + log(alpha)[q] )
      
      tau.m = NormalizeTauBrob(tsum1 )
      tau.m[is.nan(tau.m)] = 1/Q
      taus._m[[t]] = tau.m
      
      if(t>1) {
        diffs[[t]] = sum(abs( tau.m -taus._m[[t-1]]))
        
        if(   diffs[[t]]  > 100 ) break  
        if(   diffs[[t]]  < TOL  ) break
        if(t>3){
          if(  abs( (diffs[[t]] - diffs[[t-1]]  )/ diffs[[t-1]] )   < TOL ) break
        }
        print (paste('tau diff - micro iter', diffs[[t]]))
      }
    }
    
    tau[m,]  = tau_prev[m,] * (1-delta_t) + tau.m * delta_t
    alpha = colMeans(tau)
    
    TX = makeTX(tau, A1);  TY = makeTX(tau, A2)
    T0 = makeT0(tau )
    TX_Q = makeTX_PQ (tau, A1, P.00 )
    TY_Q = makeTX_PQ (tau, A2, P.00)
    
    T0_Q = makeT0_Q(tau,  P.00)
    
    xi1 =    nan_omit(TX / T0)   ;xi2 =    nan_omit(TY / T0)    ; xis = c(xi1,xi2)
    xiB1 =     nan_omit( TX_Q /  T0_Q ); xiB2 =     nan_omit( TY_Q /  T0_Q );  xiBs = c(xiB1,xiB2)
    Mu1 = MStep_Mu_Q(A1,tau) * P.11 + (1-P.11)* xi1
    Mu2 = MStep_Mu_Q(A2,tau  )* P.11+ (1-P.11)* xi2
    
    TX2 =  makeTX2(A1,tau,xi1);TY2 =  makeTX2(A2,tau,xi2)
    TX2Q = makeTX2_PQ  (tau,A1, P00 = P.00   ,xi1);  TY2Q = makeTX2_PQ  (tau,A2, P00 = P.00   ,xi2)
    
    Sigmaxi1  =   nan_omit(TX2 /T0) ;        Sigmaxi2  =   nan_omit(TY2 /T0)   
    sigmaXis = c(sqrt(Sigmaxi1),sqrt(Sigmaxi2) )
    
    SigmaxiB1  =     nan_omit(TX2Q /  T0_Q ) ;         SigmaxiB2  =     nan_omit(TY2Q /  T0_Q )
    sigmaXiBs = c(sqrt(SigmaxiB1),sqrt(SigmaxiB2) )
    
    SigmaX = MStep_SigmaX ( A1,tau,Mu1, xi1 ) * P.11 + Sigmaxi1*P.00
    SigmaY = MStep_SigmaX ( A2,tau,Mu2, xi2)* P.11 + Sigmaxi2*P.00
    
    SigmaXY = MStep_SigmaXY  (A1,A2,tau,Mu1,Mu2 ,xi1,xi2 )*P.11
    SigmaX[P.11==0] = Sigmaxi1
    SigmaY[P.11==0] = Sigmaxi2
    rhos12 =  ( nan_omit( SigmaXY  / sqrt(SigmaX * SigmaY)))
    CrossList = make_variances_Biv  (SigmaX, SigmaY,rhos12)
    
    SigmaXY = sapply( CrossList, function(x)   x[1,2])
    
    Alpha_list[[h]] = alpha
    bigtau_list[[h]] = tau
    P_list[[h]] = P.
    
    rhos = rhos12
    
    if(h>1){
      print (list('all tau', sum(abs( bigtau_list[[h]] -  bigtau_list[[h-1]]   )) ,
                  tautable(tau),P.,round(rhos,2) ))
      if (sum(abs( bigtau_list[[h]] -  bigtau_list[[h-1]] )) < TOL ) break
    }else{
      print(paste('end first round', tautable(tau)))
    }
  } 
  
  compu_time = proc.time() - tt
  
  RunList  = list(   tau=tau, Mu1=Mu1, Mu2=Mu2, P.11=P.11, rhos=rhos ,
                     SigmaXY=SigmaXY,
                     SigmaX =SigmaX, SigmaY=SigmaY,
                     xis = xis, sigmaXis=sigmaXis,
                     Alpha_list=Alpha_list, bigtau_list=bigtau_list ,P_list=P_list,
                     compu_time=compu_time)
  
  
  return(RunList)
  
}




