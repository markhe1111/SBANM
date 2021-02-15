

rm(list = ls())
  
 
library(rlist)
library(MASS)
library(Brobdingnag)
library(fpc)
  
source('~/Dropbox/DifferentialBlockModel/multiSBM/util/Tri_Utils.R')

setwd('/Users/markhe 1/Dropbox/DifferentialBlockModel/multiSBM/')
Rcpp::sourceCpp('code_Run/cpp/ClosedP_q/Bivariate/E_Biv.cpp')
Rcpp::sourceCpp('code_Run/cpp/ClosedP_q/M_Step.cpp')
Rcpp::sourceCpp('code_Run/cpp/ClosedP_q/M_Step_b.cpp')

library(randnet)
make_variances_Biv  <-function(varX_qq, varY_qq,rho_qq){
  variances = list()
  for(q in 1:Q){
    varmat = matrix(ncol=2, nrow=2, 0)
    
    varmat[1,1] = varX_qq[q]
    varmat[2,2] = varY_qq[q]
    
    varmat[1,2] = sqrt(varX_qq[q] * varY_qq[q]) * rho_qq[q]
    varmat[2,1] = sqrt(varX_qq[q] * varY_qq[q]) * rho_qq[q]
    variances[[q]] = varmat
  }
  return(variances)
}
tautable_fill  = function(tau,Q){
  
  memes = apply(tau, 1, which.max)
  q_valid  = as.numeric( names( table(memes) ) )
  qrep = rep(0,Q)
  
  for(q in 1:Q){
    if(q %in% q_valid){
      qrep[q] =     sum(memes==q)
    }else{
      qrep[q] = 0
    }
  }
  qrep
  
} 
nan_omit = function(D){
  D[is.nan(D)] = 0
  D[is.na(D)] = 0
  
  D
}
ones = function(i,Q){
  empty_vec = rep(0,Q)
  empty_vec[i] = 1
  empty_vec
}
rowMax = function(M){
  apply(M, 1, max)
}
get_EN_each_q <-function(q, A1 ,   A2 , tau ,   Mu1,  Mu2,Xis,
                         SigmaX, SigmaY, SigmaXY, sigmaXis){
  
  CovMat_q = makeCovSig( q-1, SigmaX,SigmaY , SigmaXY)  
  covNoi = diag( sigmaXis^2) 
  Mus_q = cbind(Mu1,Mu2)[q,]
  Freq_q = E_Nq_Simp_Biv   ( A1,A2, tau[,(q)],
                             Mus_q, CovMat_q,  ginv(CovMat_q),  
                             Xis,   covNoi,    ginv(covNoi))
  Freq_q
}








for(TRIAL in 1:400){
  
  Q = sample(c(3,4,5),1)
  
  meanX_qq = rnorm(Q,0,5)
  varX_qq = abs( rnorm(Q,0,5) )
  
  meanY_qq =  rnorm(Q,2,5)
  varY_qq = abs( rnorm(Q,0,5) )
  
  meanZ_qq =  rnorm(Q,4,5)
  varZ_qq = abs(rnorm(Q,0,2))
  
  corrs =   runif(Q-1,min=.01, max = .99) 
  rho_qq = c(0 , corrs)
  P.mem.fixed =c( 0  , rep(1, Q-1))
  
  # make outer products of the rhos/sigmas
  
  Xis = c(-1, 0 , 1)
  SigmaXis = abs( rnorm(Q,0,5) )
  
  meanX_qq [1] = Xis[1]
  meanY_qq [1] = Xis[2]
  meanZ_qq [1] = Xis[3]
  
  varX_qq [1] =SigmaXis[1]
  varY_qq [1] =SigmaXis[2]
  varZ_qq [1] =SigmaXis[3]
  
  simndir = rdirichlet(1,rep(1,Q))
  mixprob = .5*(simndir + rep(1/Q,Q))
  group_sizes = rmultinom(1,size = n ,mixprob)   
  #group_sizes[group_sizes==0] = 2 
  
  memberships = unlist(lapply(1:Q, function(q) rep(q, times=group_sizes[q])))
  real = cbind(P.mem.fixed , 
               meanX_qq, meanY_qq,  meanZ_qq,
               varX_qq, varY_qq, varZ_qq,
               group_sizes, rho_qq)
  
  
  variances <- make_variances_Trivariate  (varX_qq, varY_qq,varZ_qq, rho_qq)
  
  G1 = matrix(nrow=n,ncol=n,0)  
  G2 = matrix(nrow=n,ncol=n,0)
  G3 = matrix(nrow=n,ncol=n,0)
  
  for(q  in 1:Q){ 
    
    meanXq = meanX_qq[q]
    meanYq = meanY_qq[q]
    meanZq = meanZ_qq[q]
    
    var_qq = variances[[q]]                      # correlation is included in here 
    means_q = c( meanX_qq[q],  meanY_qq[q] ,  meanZ_qq[q] )
    
    mem_q = which(memberships==q)
    Not_q = setdiff(1:n,mem_q )
    
    nq = length(mem_q)
    
    nn2 = nq *( nq-1)/2
    Sim = mvtnorm::rmvnorm(n = nn2  , mean =  means_q  , sigma = var_qq    )      
    
    S1 = Sim[,1]
    S2 = Sim[,2]
    S3 = Sim[,3]
    
    b1 = fillin( nq, S1)
    b2 = fillin( nq, S2)
    b3 = fillin( nq, S3)
    
    G1[mem_q, mem_q] =  b1
    G2[mem_q, mem_q] =  b2
    G3[mem_q, mem_q] =  b3
  }
  
  # NOISE
  for(i in 1:n){
    for(j in 1:n){
      
      if(j <i){
        
        q_i = memberships[i]
        l_j = memberships[j]
        
        if(q_i != l_j   ){
          
          G1[i,j] = G1[j,i] = rnorm(n=1, mean = Xis[1], sd =  sqrt(SigmaXis[1] ))
          G2[i,j] = G2[j,i] = rnorm(n=1, mean = Xis[2], sd =  sqrt(SigmaXis[2] ))
          G3[i,j] = G3[j,i] = rnorm(n=1, mean = Xis[3], sd =  sqrt(SigmaXis[3] ))
        }
      }
    }
  } 
  
  diag(G1) = 0
  diag(G2) = 0
  diag(G3) = 0
  
  real_cors = lapply (1:Q, function(q)  cor.checker(q,memberships, G1, G2, G3))
  
  filename = paste('ReSIM_', TRIAL, '.Rdata', sep = '' )
  
  save.image( filename)
}



setwd('/Users/markhe 1/Dropbox/DifferentialBlockModel/multiSBM/data/NewSimulations/manybivSmallSim/')
files_list =  list.files()
 
file = list.files()[1]


load_and_run <-function(file){
  
  load(file)
  
  A1 = G1; A2 = G2;  
  n = ncol(A1)
  
  real = cbind(P.mem.fixed ,   meanX_qq, meanY_qq,  
               varX_qq, varY_qq,    group_sizes, rho_qq)
     
  # pretend to know the real Q
  
  runSim = function(Q, A1,  A2){
    
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
      
      initial= sapply(sphericalspectra_tot$cluster, function(x)  ones (x,Q) )
      tau = t(initial  )
      #initial = sapply(1:Q, function(x) runif(n, min=0, max=10) ) 
      #  tau = initial/rowSums(initial)
      
      P. = cbind(rep( (Q-1)/Q  ,Q), rep(1/Q,Q))
      
      #  P.[1,] = c(0,1);    
      P.11 = P.[,1];  P.00 = P.[,2]
      alpha = colMeans(tau)
      #change_P_if_noise_blocks  
      
      TX = makeTX(tau, A1)
      TY = makeTX(tau, A2)
      
      T0 = makeT0(tau )
      TX_Q = makeTX_PQ (tau, A1, P.00 )
      TY_Q = makeTX_PQ (tau, A2, P.00)
      
      T0_Q = makeT0_Q(tau,  P.00) 
      
      xi1 =   (Psi)* nan_omit(TX / T0)    + (1-Psi)* nan_omit( TX_Q /  T0_Q )
      xi2 =   (Psi)* nan_omit(TY / T0)    + (1-Psi)* nan_omit( TY_Q /  T0_Q )
      
      Xis = c(xi1,xi2)
      
     Mu1 = MStep_Mu_Q(A1,tau) * P.11 + (1-P.11)* xi1
      Mu2 = MStep_Mu_Q(A2,tau  )* P.11+ (1-P.11)* xi2
      
      
      TX2 =  makeTX2(A1,tau,xi1)
      TY2 =  makeTX2(A2,tau,xi2)
      
      # PQ is in alt
      
      TX2Q = makeTX2_PQ  (tau,A1, P00 = P.00   ,xi1)
      TY2Q = makeTX2_PQ  (tau,A2, P00 = P.00   ,xi2)
      
      Sigmaxi1  =    (Psi)* nan_omit(TX2 /T0)  +  (1-Psi)* nan_omit(TX2Q /  T0_Q )
      Sigmaxi2  =   (Psi)*  nan_omit(TY2 /T0)  +  (1-Psi)*nan_omit(TY2Q /  T0_Q )
      sigmaXis = c(sqrt(Sigmaxi1),sqrt(Sigmaxi2) )
      
      #signal part
      # NO NOISE!!!
      SigmaX = MStep_SigmaX ( A1,tau,Mu1, xi1 ) * P.11 + Sigmaxi1*P.00
      SigmaY = MStep_SigmaX ( A2,tau,Mu2, xi2)* P.11 + Sigmaxi2*P.00
      
      SigmaXY = MStep_SigmaXY  (A1,A2,tau,Mu1,Mu2 ,xi1,xi2 )*P.11
      
      rhos12 = SigmaXY  / sqrt(SigmaX * SigmaY) 
      rhos12[is.na(rhos12)] = 0 
    }
    
    TOL = 1e-300
    
    tt = proc.time()
    for(h in  1  :20 ){
      
      tau_prev = tau
      delta_t =  1- 1/(h+1) 
      #  StochSize = round( n *  delta_t^1.2 ) 
      StochSize =   min(n, round(  100+ n *  delta_t^2)  )
      #StochSize
      m= sort(sample(n,  StochSize))
      tau.m = tau[m,]
      
      QQs = sapply( 1:Q,   function(q) get_EN_each_q(q, A1 ,   A2 , tau ,  
                                                     Mu1,  Mu2, SigmaX, SigmaY, 
                                                     SigmaXY,  sigmaXis))
      
      Nq_inv1 =  as.brob( exp( QQs + log( (1-Psi)/  Psi)) )
      # Nq_inv1 =  as.brob( exp( QQs + colSums(tau) * log( (1-Psi)/  Psi)   ) )
      Nq_inv2 =   ( 1 + 1/ (Nq_inv1) )
      
      N_n = 1/ Nq_inv2
      norm_N =  nan_omit ( as.numeric( N_n/ sum(N_n) ))
      # nan_omit(norm_N)
      if(sum(norm_N)==0){
        norm_N = rep(1/Q,Q)
      }
      
      P.. =   cbind(1- norm_N ,  norm_N)
      P. = P.. * delta_t + P. * (1-delta_t)
      P.11 = P.[,1] ; P.00 = P.[,2]
      
      # Regular EM Step, w subsamp
      taus._m = list()
      
      for(t in 1:30){
        t12 = E_tau_Biv_Mat (  A1 [m,m],   A2[m,m],   tau[m,] ,  
                               Mu1,  Mu2, SigmaX, SigmaY, 
                               SigmaXY, 
                               Xis ,   sigmaXis,     P1 =P.11)
        t12[is.nan(t12)] = -Inf
        
        tsum1 = sapply(1:Q, function(q) t12[,q] + log(alpha)[q]  ) 
        #  tsum1 = sapply(1:Q, function(q) t12[,q] + log(alpha)[q]+ P.11[q]*log(Psi) +(1-P.11[q])*log(1-Psi)) 
        
        
        tau.m = NormalizeTauBrob(tsum1 )
        tau.m[is.nan(tau.m)] = 1/Q 
        taus._m[[t]] = tau.m
        
        if(t>1) {
          diffs[[t]] = sum(abs( tau.m -taus._m[[t-1]]))
          
          if(   diffs[[t]]  > 100 ) break  # this is 200 (sampled number) divided by 2
          if(   diffs[[t]]  < TOL  ) break
          if(t>3){
            if(  abs( (diffs[[t]] - diffs[[t-1]]  )/ diffs[[t-1]] )   <1e-3 ) break
          }
          print (paste('tau diff - micro iter', diffs[[t]]))
        }
      }
      
      tau[m,]  = tau_prev[m,] * (1-delta_t) + tau.m * delta_t
      alpha = colMeans(tau)
      
      #change_P_if_noise_blocks  
      
      TX = makeTX(tau, A1)
      TY = makeTX(tau, A2)
      
      T0 = makeT0(tau )
      TX_Q = makeTX_PQ (tau, A1, P.00 )
      TY_Q = makeTX_PQ (tau, A2, P.00)
      
      T0_Q = makeT0_Q(tau,  P.00) 
      
      xi1 =   (Psi)* nan_omit(TX / T0)    + (1-Psi)* nan_omit( TX_Q /  T0_Q )
      xi2 =   (Psi)* nan_omit(TY / T0)    + (1-Psi)* nan_omit( TY_Q /  T0_Q )
      
      Xis = c(xi1,xi2 )
      
      Mu1 = MStep_Mu_Q(A1,tau) * P.11 +(1-P.11)* xi1
      Mu2 = MStep_Mu_Q(A2,tau )* P.11 + (1-P.11)* xi2
      
      Mu1[P.11==0] = xi1
      Mu2[P.11==0] = xi2
      TX2 =  makeTX2(A1,tau,xi1)
      TY2 =  makeTX2(A2,tau,xi2)
      TX2Q = makeTX2_PQ  (tau,A1,P.00, xi1)
      TY2Q = makeTX2_PQ  (tau,A2,P.00, xi2)
      Sigmaxi1  =   (Psi)* nan_omit(TX2 /T0)  +   (1-Psi)* nan_omit(  (TX2Q) /  T0_Q )
      Sigmaxi2  =   (Psi)*  nan_omit(TY2 /T0)  +  (1-Psi)*nan_omit(   (TY2Q) /  T0_Q )
      
      sigmaXis = c(sqrt(Sigmaxi1),sqrt(Sigmaxi2) )
      
      #signal part
      # NO NOISE!!!
      SigmaX = MStep_SigmaX ( A1,tau,Mu1, xi1 )*P.11 + Sigmaxi1*P.00
      SigmaY = MStep_SigmaX ( A2,tau,Mu2, xi2)* P.11 + Sigmaxi2*P.00
      
      SigmaXY = MStep_SigmaXY  (A1,A2,tau,Mu1,Mu2 ,xi1,xi2 )*P.11
      
      SigmaX[P.11==0] = Sigmaxi1
      SigmaY[P.11==0] = Sigmaxi2
      rhos12 =  ( nan_omit( SigmaXY  / sqrt(SigmaX * SigmaY))) 
      
      rhosmax =  rhos12
      rhosmax[rhosmax<0] = 0
      
      
      # make this rowmax or rowmeans
      # then make into zero if 
      CrossList = make_variances_Biv  (SigmaX, SigmaY,rhosmax)
      
      SigmaXY = sapply( CrossList, function(x)   x[1,2])
      
      Alpha_list[[h]] = alpha
      bigtau_list[[h]] = tau
      P_list[[h]] = P.
      
      #ELBOs[[h]] =  try(trivELBO (  A1,   A2,  A3, tau,     Mu1,  Mu2,  Mu3, SigmaX, SigmaY, SigmaZ,
      #                              SigmaXY, SigmaYZ,  SigmaXZ,  Xis ,   sigmaXis,   P1 =P.11))
      
      rhos = rhosmax
      
      if(h>1){
        print (list('all tau', sum(abs( bigtau_list[[h]] -  bigtau_list[[h-1]]   )) , 
                    tautable(tau),P.,round(rhos,2) ))
        if (sum(abs( bigtau_list[[h]] -  bigtau_list[[h-1]] )) < 1e-3 ) break
      }else{
        print(paste('end first round', tautable(tau)))
      }
    }  
    compu_time = proc.time() - tt
    
    RunList  = list(   tau=tau, Mu1=Mu1, Mu2=Mu2, P.11=P.11, rhos=rhos ,
                       SigmaXY=SigmaXY, 
                       SigmaX =SigmaX, SigmaY=SigmaY, 
                       Xis = Xis, sigmaXis=sigmaXis,
                       Alpha_list=Alpha_list, bigtau_list=bigtau_list ,P_list=P_list, ELBOS=ELBOs, compu_time=compu_time)
    
    
    return(RunList)
  }
  result = runSim(Q,A1,  A2 )
  result
}

Run_1_10 = lapply(1:10 , function(XX)  try( load_and_run (files_list[XX]  )  ))

compu_tot2 = proc.time()
for(R in 1:50 ){
  Runs [[R]]  = try( load_and_run (files_list[R] )) 
}
comptime_2 = proc.time() - compu_tot2

 

setwd('~/Dropbox/DifferentialBlockModel/multiSBM/Results/Sim_results/current/ManySmallBiv')
#save.image('manysmall_1_50.Rdata')
#load('manysmall_1_50.Rdata')




#setwd('~/Dropbox/DifferentialBlockModel/multiSBM/Da')
setwd("~/Dropbox/DifferentialBlockModel/multiSBM/data/NewSimulations/manybivSmallSim/")
files_list = list.files()

symdiff <- function (s1, s2) {
  return(union(setdiff(s1, s2), setdiff(s2, s1)))
}
jaccard <- function (s1, s2) {
  return(length(symdiff(s1, s2)) / length(union(s1, s2)))
}

R_Agreememt = list()

for(R in 1:50 ){
 
  RunR = Runs[[R]]
  load(files_list[R])
  
  if( class(RunR) == "try-error" ) {
    R_Agreememt[[R]] <- 'error'
  }else{
    sorted_tautable = sort ( tautable( RunR$tau))
    sorted_realtable = sort( table(memberships) )
    R_Agreememt[[R]] = jaccard (sorted_tautable, sorted_realtable)
  }
}



unlist( R_Agreememt )

reRUn_List = which ( sapply (R_Agreememt   ,  function(x) x !=0   ) )


 

re_RUN <-function(file){
  
  load(file)
  
  A1 = G1; A2 = G2;  
  n = ncol(A1)
  
  real = cbind(P.mem.fixed ,   meanX_qq, meanY_qq,  
               varX_qq, varY_qq,    group_sizes, rho_qq)
  
  # pretend to know the real Q
  
  runSim = function(Q, A1,  A2){
    
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
      
      initial= sapply(sphericalspectra_tot$cluster, function(x)  ones (x,Q) )
      tau = t(initial  )
      #initial = sapply(1:Q, function(x) runif(n, min=0, max=10) ) 
      #  tau = initial/rowSums(initial)
      
      P. = cbind(rep( (Q-1)/Q  ,Q), rep(1/Q,Q))
      
      #  P.[1,] = c(0,1);    
      P.11 = P.[,1];  P.00 = P.[,2]
      alpha = colMeans(tau)
      #change_P_if_noise_blocks  
      
      TX = makeTX(tau, A1)
      TY = makeTX(tau, A2)
      
      T0 = makeT0(tau )
      TX_Q = makeTX_PQ (tau, A1, P.00 )
      TY_Q = makeTX_PQ (tau, A2, P.00)
      
      T0_Q = makeT0_Q(tau,  P.00) 
      
      xi1 =   (Psi)* nan_omit(TX / T0)    + (1-Psi)* nan_omit( TX_Q /  T0_Q )
      xi2 =   (Psi)* nan_omit(TY / T0)    + (1-Psi)* nan_omit( TY_Q /  T0_Q )
      
      Xis = c(xi1,xi2)
      
      Mu1 = MStep_Mu_Q(A1,tau) * P.11 + (1-P.11)* xi1
      Mu2 = MStep_Mu_Q(A2,tau  )* P.11+ (1-P.11)* xi2
      
      
      TX2 =  makeTX2(A1,tau,xi1)
      TY2 =  makeTX2(A2,tau,xi2)
      
      # PQ is in alt
      
      TX2Q = makeTX2_PQ  (tau,A1, P00 = P.00   ,xi1)
      TY2Q = makeTX2_PQ  (tau,A2, P00 = P.00   ,xi2)
      
      Sigmaxi1  =    (Psi)* nan_omit(TX2 /T0)  +  (1-Psi)* nan_omit(TX2Q /  T0_Q )
      Sigmaxi2  =   (Psi)*  nan_omit(TY2 /T0)  +  (1-Psi)*nan_omit(TY2Q /  T0_Q )
      sigmaXis = c(sqrt(Sigmaxi1),sqrt(Sigmaxi2) )
      
      #signal part
      # NO NOISE!!!
      SigmaX = MStep_SigmaX ( A1,tau,Mu1, xi1 ) * P.11 + Sigmaxi1*P.00
      SigmaY = MStep_SigmaX ( A2,tau,Mu2, xi2)* P.11 + Sigmaxi2*P.00
      
      SigmaXY = MStep_SigmaXY  (A1,A2,tau,Mu1,Mu2 ,xi1,xi2 )*P.11
      
      rhos12 = SigmaXY  / sqrt(SigmaX * SigmaY) 
      rhos12[is.na(rhos12)] = 0 
    }
    
    TOL = 1e-30
    
    tt = proc.time()
    for(h in  1  :40 ){
      
      tau_prev = tau
      delta_t =  1- 1/(h+1) 
      #  StochSize = round( n *  delta_t^1.2 ) 
      StochSize =   min(n, round(  100+ n *  delta_t^2)  )
      #StochSize
      m= sort(sample(n,  StochSize))
      tau.m = tau[m,]
      
      QQs = sapply( 1:Q,   function(q) get_EN_each_q(q, A1 ,   A2 , tau ,   
                                                     Mu1,  Mu2, Xis,
                                                     SigmaX, SigmaY, 
                                                     SigmaXY,  sigmaXis))
      
      Nq_inv1 =  as.brob( exp( QQs + log( (1-Psi)/  Psi)) )
      # Nq_inv1 =  as.brob( exp( QQs + colSums(tau) * log( (1-Psi)/  Psi)   ) )
      Nq_inv2 =   ( 1 + 1/ (Nq_inv1) )
      
      N_n = 1/ Nq_inv2
      norm_N =  nan_omit ( as.numeric( N_n/ sum(N_n) ))
      # nan_omit(norm_N)
      if(sum(norm_N)==0){
        norm_N = rep(1/Q,Q)
      }
      
      P.. =   cbind(1- norm_N ,  norm_N)
      P. = P.. * delta_t + P. * (1-delta_t)
      P.11 = P.[,1] ; P.00 = P.[,2]
      
      # Regular EM Step, w subsamp
      taus._m = list()
      
      for(t in 1:30){
        t12 = E_tau_Biv_Mat (  A1 [m,m],   A2[m,m],   tau[m,] ,  
                               Mu1,  Mu2, SigmaX, SigmaY, 
                               SigmaXY, 
                               Xis ,   sigmaXis,     P1 =P.11)
        t12[is.nan(t12)] = -Inf
        
        tsum1 = sapply(1:Q, function(q) t12[,q] + log(alpha)[q]  ) 
        #  tsum1 = sapply(1:Q, function(q) t12[,q] + log(alpha)[q]+ P.11[q]*log(Psi) +(1-P.11[q])*log(1-Psi)) 
        
        
        tau.m = NormalizeTauBrob(tsum1 )
        tau.m[is.nan(tau.m)] = 1/Q 
        taus._m[[t]] = tau.m
        
        if(t>1) {
          diffs[[t]] = sum(abs( tau.m -taus._m[[t-1]]))
          
          if(   diffs[[t]]  > 100 ) break  # this is 200 (sampled number) divided by 2
          if(   diffs[[t]]  < TOL  ) break
          if(t>3){
            if(  abs( (diffs[[t]] - diffs[[t-1]]  )/ diffs[[t-1]] )   <1e-3 ) break
          }
          print (paste('tau diff - micro iter', diffs[[t]]))
        }
      }
      
      tau[m,]  = tau_prev[m,] * (1-delta_t) + tau.m * delta_t
      alpha = colMeans(tau)
      
      #change_P_if_noise_blocks  
      
      TX = makeTX(tau, A1)
      TY = makeTX(tau, A2)
      
      T0 = makeT0(tau )
      TX_Q = makeTX_PQ (tau, A1, P.00 )
      TY_Q = makeTX_PQ (tau, A2, P.00)
      
      T0_Q = makeT0_Q(tau,  P.00) 
      
      xi1 =   (Psi)* nan_omit(TX / T0)    + (1-Psi)* nan_omit( TX_Q /  T0_Q )
      xi2 =   (Psi)* nan_omit(TY / T0)    + (1-Psi)* nan_omit( TY_Q /  T0_Q )
      
      Xis = c(xi1,xi2 )
      
      Mu1 = MStep_Mu_Q(A1,tau) * P.11 +(1-P.11)* xi1
      Mu2 = MStep_Mu_Q(A2,tau )* P.11 + (1-P.11)* xi2
      
      Mu1[P.11==0] = xi1
      Mu2[P.11==0] = xi2
      TX2 =  makeTX2(A1,tau,xi1)
      TY2 =  makeTX2(A2,tau,xi2)
      TX2Q = makeTX2_PQ  (tau,A1,P.00, xi1)
      TY2Q = makeTX2_PQ  (tau,A2,P.00, xi2)
      Sigmaxi1  =   (Psi)* nan_omit(TX2 /T0)  +   (1-Psi)* nan_omit(  (TX2Q) /  T0_Q )
      Sigmaxi2  =   (Psi)*  nan_omit(TY2 /T0)  +  (1-Psi)*nan_omit(   (TY2Q) /  T0_Q )
      
      sigmaXis = c(sqrt(Sigmaxi1),sqrt(Sigmaxi2) )
      
      #signal part
      # NO NOISE!!!
      SigmaX = MStep_SigmaX ( A1,tau,Mu1, xi1 )*P.11 + Sigmaxi1*P.00
      SigmaY = MStep_SigmaX ( A2,tau,Mu2, xi2)* P.11 + Sigmaxi2*P.00
      
      SigmaXY = MStep_SigmaXY  (A1,A2,tau,Mu1,Mu2 ,xi1,xi2 )*P.11
      
      SigmaX[P.11==0] = Sigmaxi1
      SigmaY[P.11==0] = Sigmaxi2
      rhos12 =  ( nan_omit( SigmaXY  / sqrt(SigmaX * SigmaY))) 
      
      rhosmax =  rhos12
      rhosmax[rhosmax<0] = 0
      
      
      # make this rowmax or rowmeans
      # then make into zero if 
      CrossList = make_variances_Biv  (SigmaX, SigmaY,rhosmax)
      
      SigmaXY = sapply( CrossList, function(x)   x[1,2])
      
      Alpha_list[[h]] = alpha
      bigtau_list[[h]] = tau
      P_list[[h]] = P.
       
      rhos = rhosmax
      
      if(h>1){
        print (list('all tau', sum(abs( bigtau_list[[h]] -  bigtau_list[[h-1]]   )) , 
                    tautable(tau),P.,round(rhos,2) ))
        if (sum(abs( bigtau_list[[h]] -  bigtau_list[[h-1]] )) < 1e-3 ) break
      }else{
        print(paste('end first round', tautable(tau)))
      }
    }  
    compu_time = proc.time() - tt
    
    RunList  = list(   tau=tau, Mu1=Mu1, Mu2=Mu2, P.11=P.11, rhos=rhos ,
                       SigmaXY=SigmaXY, 
                       SigmaX =SigmaX, SigmaY=SigmaY, 
                       Xis = Xis, sigmaXis=sigmaXis,
                       Alpha_list=Alpha_list, bigtau_list=bigtau_list ,P_list=P_list, ELBOS=ELBOs, compu_time=compu_time)
    
    
    return(RunList)
  }
  result = runSim(Q,A1,  A2 )
  result
}

 




compu_totRe = proc.time()
for(R in reRUn2 ){
  Runs [[R]]  = try( re_RUN (files_list[R] )) 
}
ReCompu_totRe = proc.time() - compu_totRe

 

Agree2 = list()
for(R in 1:50 ){
  RunR = Runs[[R]]
  load(files_list[R])
  
  if( class(RunR) == "try-error" ) {
    Agree2[[R]] <- 'error'
  }else{
    sorted_tautable = sort ( tautable( RunR$tau))
    sorted_realtable = sort( table(memberships) )
    Agree2[[R]] = jaccard (sorted_tautable, sorted_realtable)
  }
}
unlist( Agree2 )
reRUn2 = which ( sapply (Agree2   ,  function(x) x !=0   ) )
length( reRUn2 )
setwd('~/Dropbox/DifferentialBlockModel/multiSBM/Results/Sim_results/current/ManySmallBiv/')
load('~/Dropbox/DifferentialBlockModel/multiSBM/Results/Sim_results/current/ManySmallBiv/Reran_31_Worked.Rdata')


setwd('/Users/markhe 1/Dropbox/DifferentialBlockModel/multiSBM/')
Rcpp::sourceCpp('code_Run/cpp/ClosedP_q/Bivariate/E_Biv.cpp')
Rcpp::sourceCpp('code_Run/cpp/ClosedP_q/M_Step.cpp')
Rcpp::sourceCpp('code_Run/cpp/ClosedP_q/M_Step_b.cpp')

library(randnet)
library(MASS)
#setwd('~/Dropbox/DifferentialBlockModel/multiSBM/Da')
setwd("~/Dropbox/DifferentialBlockModel/multiSBM/data/NewSimulations/manybivSmallSim/")
files_list = list.files()
compu_totRe = proc.time()
for(R in reRUn2 ){
  Runs [[R]]  = try( re_RUN (files_list[R] )) 
}
ReCompu_totRe = proc.time() - compu_totRe



Agree3 = list()
for(R in 1:50 ){
  RunR = Runs[[R]]
  load(files_list[R])
  
  if( class(RunR) == "try-error" ) {
    Agree3[[R]] <- 'error'
  }else{
    sorted_tautable = sort ( tautable( RunR$tau))
    sorted_realtable = sort( table(memberships) )
    Agree3[[R]] = jaccard (sorted_tautable, sorted_realtable)
  }
}

unlist( Agree3 )
reRUn3 = which ( sapply (Agree3   ,  function(x) x !=0   ) )
length( reRUn3 )

load('~/Dropbox/DifferentialBlockModel/multiSBM/Results/Sim_results/current/ManySmallBiv/Reran_35_Worked.Rdata')


re_RUN3 <-function(file){
  
  load(file)
  
  A1 = G1; A2 = G2;  
  n = ncol(A1)
  
  real = cbind(P.mem.fixed ,   meanX_qq, meanY_qq,  
               varX_qq, varY_qq,    group_sizes, rho_qq)
  
  # pretend to know the real Q
  
  runSim = function(Q, A1,  A2){
    
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
      
      initial= sapply(sphericalspectra_tot$cluster, function(x)  ones (x,Q) )
      tau = t(initial  )
      #initial = sapply(1:Q, function(x) runif(n, min=0, max=10) ) 
      #  tau = initial/rowSums(initial)
      
      P. = cbind(rep( (Q-1)/Q  ,Q), rep(1/Q,Q))
      
      #  P.[1,] = c(0,1);    
      P.11 = P.[,1];  P.00 = P.[,2]
      alpha = colMeans(tau)
      #change_P_if_noise_blocks  
      
      TX = makeTX(tau, A1)
      TY = makeTX(tau, A2)
      
      T0 = makeT0(tau )
      TX_Q = makeTX_PQ (tau, A1, P.00 )
      TY_Q = makeTX_PQ (tau, A2, P.00)
      
      T0_Q = makeT0_Q(tau,  P.00) 
      
      xi1 =   (Psi)* nan_omit(TX / T0)    + (1-Psi)* nan_omit( TX_Q /  T0_Q )
      xi2 =   (Psi)* nan_omit(TY / T0)    + (1-Psi)* nan_omit( TY_Q /  T0_Q )
      
      Xis = c(xi1,xi2)
      
      Mu1 = MStep_Mu_Q(A1,tau) * P.11 + (1-P.11)* xi1
      Mu2 = MStep_Mu_Q(A2,tau  )* P.11+ (1-P.11)* xi2
      
      TX2 =  makeTX2(A1,tau,xi1)
      TY2 =  makeTX2(A2,tau,xi2)
      
      # PQ is in alt
      
      TX2Q = makeTX2_PQ  (tau,A1, P00 = P.00   ,xi1)
      TY2Q = makeTX2_PQ  (tau,A2, P00 = P.00   ,xi2)
      
      Sigmaxi1  =    (Psi)* nan_omit(TX2 /T0)  +  (1-Psi)* nan_omit(TX2Q /  T0_Q )
      Sigmaxi2  =   (Psi)*  nan_omit(TY2 /T0)  +  (1-Psi)*nan_omit(TY2Q /  T0_Q )
      sigmaXis = c(sqrt(Sigmaxi1),sqrt(Sigmaxi2) )
      
      #signal part
      # NO NOISE!!!
      SigmaX = MStep_SigmaX ( A1,tau,Mu1, xi1 ) * P.11 + Sigmaxi1*P.00
      SigmaY = MStep_SigmaX ( A2,tau,Mu2, xi2)* P.11 + Sigmaxi2*P.00
      
      SigmaXY = MStep_SigmaXY  (A1,A2,tau,Mu1,Mu2 ,xi1,xi2 )*P.11
      
      rhos12 = SigmaXY  / sqrt(SigmaX * SigmaY) 
      rhos12[is.na(rhos12)] = 0 
    }
    
    TOL = 1e-30
    
    tt = proc.time()
    for(h in  1  :40 ){
      
      tau_prev = tau
      delta_t =  1- 1/(h+1) 
      #  StochSize = round( n *  delta_t^1.2 ) 
      StochSize =   min(n, round(  300+ n *  delta_t^2)  )
      #StochSize
      m= sort(sample(n,  StochSize))
      tau.m = tau[m,]
      
      QQs = sapply( 1:Q,   function(q) get_EN_each_q(q, A1 ,   A2 , tau ,   
                                                     Mu1,  Mu2, Xis,
                                                     SigmaX, SigmaY, 
                                                     SigmaXY,  sigmaXis))
      
      Nq_inv1 =  as.brob( exp( QQs + log( (1-Psi)/  Psi)) )
      # Nq_inv1 =  as.brob( exp( QQs + colSums(tau) * log( (1-Psi)/  Psi)   ) )
      Nq_inv2 =   ( 1 + 1/ (Nq_inv1) )
      
      N_n = 1/ Nq_inv2
      norm_N =  nan_omit ( as.numeric( N_n/ sum(N_n) ))
      # nan_omit(norm_N)
      if(sum(norm_N)==0){
        norm_N = rep(1/Q,Q)
      }
      
      P.. =   cbind(1- norm_N ,  norm_N)
      P. = P.. * delta_t + P. * (1-delta_t)
      P.11 = P.[,1] ; P.00 = P.[,2]
      
      # Regular EM Step, w subsamp
      taus._m = list()
      
      for(t in 1:30){
        t12 = E_tau_Biv_Mat (  A1 [m,m],   A2[m,m],   tau[m,] ,  
                               Mu1,  Mu2, SigmaX, SigmaY, 
                               SigmaXY, 
                               Xis ,   sigmaXis,     P1 =P.11)
        t12[is.nan(t12)] = -Inf
        
        tsum1 = sapply(1:Q, function(q) t12[,q] + log(alpha)[q]  ) 
        #  tsum1 = sapply(1:Q, function(q) t12[,q] + log(alpha)[q]+ P.11[q]*log(Psi) +(1-P.11[q])*log(1-Psi)) 
        
        
        tau.m = NormalizeTauBrob(tsum1 )
        tau.m[is.nan(tau.m)] = 1/Q 
        taus._m[[t]] = tau.m
        
        if(t>1) {
          diffs[[t]] = sum(abs( tau.m -taus._m[[t-1]]))
          
          if(   diffs[[t]]  > 100 ) break  # this is 200 (sampled number) divided by 2
          if(   diffs[[t]]  < TOL  ) break
          if(t>3){
            if(  abs( (diffs[[t]] - diffs[[t-1]]  )/ diffs[[t-1]] )   <1e-3 ) break
          }
          print (paste('tau diff - micro iter', diffs[[t]]))
        }
      }
      
      tau[m,]  = tau_prev[m,] * (1-delta_t) + tau.m * delta_t
      alpha = colMeans(tau)
      
      #change_P_if_noise_blocks  
      
      TX = makeTX(tau, A1)
      TY = makeTX(tau, A2)
      
      T0 = makeT0(tau )
      TX_Q = makeTX_PQ (tau, A1, P.00 )
      TY_Q = makeTX_PQ (tau, A2, P.00)
      
      T0_Q = makeT0_Q(tau,  P.00) 
      
      xi1 =   (Psi)* nan_omit(TX / T0)    + (1-Psi)* nan_omit( TX_Q /  T0_Q )
      xi2 =   (Psi)* nan_omit(TY / T0)    + (1-Psi)* nan_omit( TY_Q /  T0_Q )
      
      Xis = c(xi1,xi2 )
      
      Mu1 = MStep_Mu_Q(A1,tau) * P.11 +(1-P.11)* xi1
      Mu2 = MStep_Mu_Q(A2,tau )* P.11 + (1-P.11)* xi2
      
      Mu1[P.11==0] = xi1
      Mu2[P.11==0] = xi2
      TX2 =  makeTX2(A1,tau,xi1)
      TY2 =  makeTX2(A2,tau,xi2)
      TX2Q = makeTX2_PQ  (tau,A1,P.00, xi1)
      TY2Q = makeTX2_PQ  (tau,A2,P.00, xi2)
      Sigmaxi1  =   (Psi)* nan_omit(TX2 /T0)  +   (1-Psi)* nan_omit(  (TX2Q) /  T0_Q )
      Sigmaxi2  =   (Psi)*  nan_omit(TY2 /T0)  +  (1-Psi)*nan_omit(   (TY2Q) /  T0_Q )
      
      sigmaXis = c(sqrt(Sigmaxi1),sqrt(Sigmaxi2) )
      
      #signal part
      # NO NOISE!!!
      SigmaX = MStep_SigmaX ( A1,tau,Mu1, xi1 )*P.11 + Sigmaxi1*P.00
      SigmaY = MStep_SigmaX ( A2,tau,Mu2, xi2)* P.11 + Sigmaxi2*P.00
      
      SigmaXY = MStep_SigmaXY  (A1,A2,tau,Mu1,Mu2 ,xi1,xi2 )*P.11
      
      SigmaX[P.11==0] = Sigmaxi1
      SigmaY[P.11==0] = Sigmaxi2
      rhos12 =  ( nan_omit( SigmaXY  / sqrt(SigmaX * SigmaY))) 
      
      rhosmax =  rhos12
      rhosmax[rhosmax<0] = 0
      
      
      # make this rowmax or rowmeans
      # then make into zero if 
      CrossList = make_variances_Biv  (SigmaX, SigmaY,rhosmax)
      
      SigmaXY = sapply( CrossList, function(x)   x[1,2])
      
      Alpha_list[[h]] = alpha
      bigtau_list[[h]] = tau
      P_list[[h]] = P.
      
      rhos = rhosmax
      
      if(h>1){
        print (list('all tau', sum(abs( bigtau_list[[h]] -  bigtau_list[[h-1]]   )) , 
                    tautable(tau),P.,round(rhos,2) ))
        if (sum(abs( bigtau_list[[h]] -  bigtau_list[[h-1]] )) < 1e-6 ) break
      }else{
        print(paste('end first round', tautable(tau)))
      }
    }  
    compu_time = proc.time() - tt
    RunList  = list(   tau=tau, Mu1=Mu1, Mu2=Mu2, P.11=P.11, rhos=rhos ,
                       SigmaXY=SigmaXY, 
                       SigmaX =SigmaX, SigmaY=SigmaY, 
                       Xis = Xis, sigmaXis=sigmaXis,
                       Alpha_list=Alpha_list, bigtau_list=bigtau_list ,P_list=P_list, ELBOS=ELBOs, compu_time=compu_time)
    
    
    return(RunList)
  }
  result = runSim(Q,A1,  A2 )
  result
}



#setwd('~/Dropbox/DifferentialBlockModel/multiSBM/Da')
setwd("~/Dropbox/DifferentialBlockModel/multiSBM/data/NewSimulations/manybivSmallSim/")
files_list = list.files()
compu_totRe = proc.time()
for(R in c(42,50) ){
  Runs [[R]]  = try( re_RUN3 (files_list[R] )) 
}
ReCompu_totRe = proc.time() - compu_totRe



 
get_agreement = function(Runs){
  
  symdiff <- function (s1, s2) {
    return(union(setdiff(s1, s2), setdiff(s2, s1)))
  }
  jaccard <- function (s1, s2) {
    return(length(symdiff(s1, s2)) / length(union(s1, s2)))
  }
  
  setwd("~/Dropbox/DifferentialBlockModel/multiSBM/data/NewSimulations/manybivSmallSim/")
  files_list = list.files()
  #list.files()
  AgreeMent = list()
  
  for(R in 1:50 ){
    RunR = Runs[[R]]
    load(files_list[R])
    
    if( class(RunR) == "list"){
      sorted_tautable = sort ( tautable( RunR$tau))
      sorted_realtable = sort( table(memberships) )
      AgreeMent[[R]] = jaccard (sorted_tautable, sorted_realtable)
    }else{
      AgreeMent[[R]] <- 'error'
    }
  }
  AgreeMent
}
Agreements = get_agreement(Runs)
sum( unlist(Agreements)==0 )

#('~/Dropbox/DifferentialBlockModel/multiSBM/Results/Sim_results/current/ManySmallBiv/Reran_31_Worked.Rdata')
load('~/Dropbox/DifferentialBlockModel/multiSBM/Results/Sim_results/current/ManySmallBiv/Reran_35_Worked.Rdata')


not_agree = which( unlist(Agreements)!=0)
agree = which( unlist(Agreements)!=0)



library(fpc)

make_cov2 <-function(SigmaX_i, SigmaY_i, SigmaXY_i  ){
  covar = matrix(ncol=2, nrow=2,0)
  covar[1,1] = SigmaX_i
  covar[1,2] = covar[2,1] = SigmaXY_i
  covar[2,2] = SigmaY_i 
  covar
}


 setwd('/Users/markhe 1/Dropbox/DifferentialBlockModel/multiSBM/data/NewSimulations/manybivSmallSim/')
files_list =  list.files() 


file = files_list[XX]

load_and_bhattachar <-function(file){
  
  load(file)
  Q 
  real
  
  est_memberships = apply(tau, 1, which.max)
  
  
  
  cov_mats = lapply(1:Q, function(q)    make_cov2 (  SigmaX[q] ,SigmaY[q],    SigmaXY[q]   )   )
   
  rhos12 = rhos 
  #Mus_real  = cbind(real,Mu2)
  
  mu_list = lapply(1:Q, function(q)  Mus[q,])
  
  Bdist_Mat = matrix(nrow = Q, ncol = Q, 0)
  colnames(Bdist_Mat) = paste('B_d', 1:Q)
  
  for(q1 in  1:Q){
    
    for(q2 in 1:Q){
      Bdist_Mat [q2,q1]=   bhattacharyya.dist(  Mus [q1,]  , Mus [q2,],
                                                Sigma1 =   cov_mats[[q1]]  , 
                                                Sigma2 =   cov_mats[[q2]] )
    }
  } 
  
  Bdist_Mat
  
}



#setwd('/Users/markhe 1/Dropbox/DifferentialBlockModel/multiSBM/data/NewSimulations/manySmallSim')
files_list =  list.files() 

Agree_Bhatta  = list()
for(XX in agree ){
  Agree_Bhatta[[XX]] =   load_and_bhattachar (files_list[XX])
} 

make_mu_table <-function(Res  ){
  
  tau = Res$tau
  Q = ncol(tau)
  Samp = Res$Samp
  
  Mu_x = Res$Mu1 
  Mu_y = Res$Mu2 
  Mu_z = Res$Mu3 
  P.11 = Res$P.11
  
  SigmaX = Res$SigmaX 
  SigmaY = Res$SigmaY 
  SigmaZ = Res$SigmaZ 
  
  SigmaXY = Res$SigmaXY
  SigmaYZ = Res$SigmaYZ 
  SigmaXZ = Res$SigmaXZ
  rhos = Res$rhos
  
  make_cov3 <-function(SigmaX_i, SigmaY_i,SigmaZ_i, SigmaXY_i,SigmaYZ_i,SigmaXZ_i ){
    covar = matrix(ncol=3, nrow=3,0)
    covar[1,1] = SigmaX_i
    covar[1,2] = covar[2,1] = SigmaXY_i
    covar[2,2] = SigmaY_i
    covar[2,3] = covar[3,2] = SigmaYZ_i
    covar[3,3] = SigmaZ_i
    covar[1,3] = covar[3,1] = SigmaXZ_i
    covar
  }
  memberships = apply(tau, 1, which.max)
  
  cov_mats = lapply(1:Q, function(q)    make_cov3 (  SigmaX[q] ,SigmaY[q], SigmaZ[q],     SigmaXY[q] ,SigmaYZ[q], SigmaXZ[q] )    )
  rhos12 = rhos[,1]
  
  Mus = cbind(Mu_x,Mu_y,Mu_z)
  mu_list = lapply(1:Q, function(q)  Mus[q,])
  
  
  
  resulting_table  =   data.frame(index =seq(Q), 
                                  SigProb = P.11,
                                  n=tautable_fill(tau, Q) ,
                                  'rho'=rhos12,  #rhos23, rhos13 ,
                                  mu_x = Mu_x ,sigma_x=sqrt(SigmaX), 
                                  mu_y = Mu_y, sigma_y=sqrt(SigmaY),
                                  mu_z = Mu_z ,sigma_z=sqrt(SigmaZ)                                  )
  
  
  sorted_table = resulting_table[ order( resulting_table$rho),] 
  
  
  
  Bdist_Mat = matrix(nrow = Q, ncol = Q, 0)
  colnames(Bdist_Mat) = paste('B_d', 1:Q)
  
  for(q1 in  sorted_table$index){
    for(q2 in  sorted_table$index){
      Bdist_Mat [q2,q1]=   bhattacharyya.dist(  Mus [q1,]  , Mus [q2,],
                                                Sigma1 =   cov_mats[[q1]]  , Sigma2 =   cov_mats[[q2]] )
    }
  } 
  ##race
  #1= White; 2= Black; 3= Native American/Alaska Native; 4=Asian; 5=More Than One Race;
  # 6=Hawaiian/Pacific; 9=Unknown/Unreported
  ## race2
  #1= White; 2= Black; 3= others;
  ##sex
  #1=male, 2=female
  ##ethnicity
  #1 = Hispanic/latino, 2= Non-hispanic-Latino\
  # table_B
  
  final_table = cbind( sorted_table, Bdist_Mat)
  
  
  table_B = round( final_table [final_table$n>0 , ] ,2)
  #table_B$mu_1 = paste(table_B$mu_1, "pm",sep='' )
  #table_B$mu_2 = paste(table_B$mu_2, "pm" ,sep='')
  #table_B$mu_3 = paste(table_B$mu_3, "pm",sep='' )
  # Tsig_order = table_B[order(table_B$SigProb ),]
  #  Trho_order = Tsig_order[order(Tsig_order$rho ),]
  #  Trho_order
  table_B
}






