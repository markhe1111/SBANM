get_mahala_ginv <-function(dist, covv ){
  
  
  invcovdist =ginv(covv) %*% t(dist)
  mahala_out = dist%*% invcovdist
  mahala_out
}
get_mahala_inv <-function(dist, covv ){
  
  invcovdist = solve(covv) %*% t(dist)
  mahala_out = dist%*% invcovdist
  mahala_out
}

   
 
make_cov2 <-function(SigmaX_i, SigmaY_i, SigmaXY_i ){
  covar = matrix(ncol=2, nrow=2,0)
  covar[1,1] = SigmaX_i
  covar[1,2] = covar[2,1] = SigmaXY_i
  covar[2,2] = SigmaY_i
  covar
}

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

bh_reject <- function (pvals, alpha) {
  pvals_adj <- length(pvals) * pvals / rank(pvals)
  
  sum.pval = sum(pvals_adj <= alpha) 
  if(is.na(sum.pval)) sum.pval <- 0
  
  if (sum.pval > 0) {
    thres <- max(pvals[pvals_adj <= alpha])
    return(which(pvals <= thres))
  } else {
    return(integer(0))
  }
}

bhy <-  function(pvals, alpha = 0.05){
    
    # Sorted p-vals
    sp = sort(pvals)
    
    # Save original order of p-vals
    ord = order(pvals)
    
    # Find bhy cutoff
    nums = 1:length(pvals)
    cms = cumsum(1/nums)
    
    # Find which p-vals are less than bh cutoff
    under = sp < (nums/(length(pvals)*cms)*alpha)
    
    # Return indices of significant p-vals
    if(sum(under) == 0){
      return(c())
    }else{
      cutoff = max(which(under))
      return(ord[1:cutoff])
    }
  }

 
tautable = function(tau){
  table(apply(tau, 1, which.max))
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
 


cor.checker <-function(q, memberships , G1=A1, G2=A2, G3=A3){
  q_1 =  G1[ memberships==q, memberships==q]
  q_2 =  G2[ memberships==q, memberships==q]
  q_3 =  G3[ memberships==q, memberships==q]
  cor12 =  cor(    q_1[ lower.tri (q_1)],  q_2[ lower.tri (q_2)])
  cor23 =  cor(    q_2[ lower.tri (q_2)],  q_3[ lower.tri (q_3)])
  cor13 =  cor(    q_1[ lower.tri (q_1)],  q_3[ lower.tri (q_3)])
  
  c(cor12, cor23, cor13)
}

NormalizeTauBrob <-function( logtau0){
  Q = ncol(logtau0)
  n = nrow(logtau0)
  ninorm = matrix(nrow=n,ncol=Q,0)
  
  for(i in 1:n){
    
    ni =list()
    for( q in 1:Q){
      ni[[q]] = exp( as.brob( (logtau0[i,q])))
    }
    nisum = 0
    for( q in 2:(Q+1)){
      nisum = nisum +ni[[q-1]]
    }
    
    
    
    for( q in 1:Q){
      ni[[q]]/nisum
      ninorm[i,q] =  as.numeric( as.brob(ni[[q]]/nisum))
    }
  }
  return(ninorm)
}
dont <-function(){
  Bdist_XY <-function (q,    Mu1,Mu2,  SigmaX,SigmaY , SigmaXY,   xi, delta,sigmaxi, sigmadelta, rho_xi) {
    dist_q = t( c( Mu1[q] - xi ,  
                   Mu2[q] - delta))
    cov_q = make_cov(SigmaX[q],SigmaY[q], SigmaXY[q])
    covv_noise =  rho_xi*sigmaxi * sigmadelta
    cov_xi = make_cov( sigmaxi^2, sigmadelta^2 ,  covv_noise )
    get_bhatta(dist_q, cov_q , cov_xi)
    
  }
  get_TwoSampleTestStat <-function(dist, SampleSigma){
    d1 = get_mahala(dist, SampleSigma)
    return(d1 )
  }
  Test_XY <-function (q,   Mu1,Mu2,  SigmaX,SigmaY , SigmaXY,    xi, delta,sigmaxi, sigmadelta, rho_xi) {
    
    dist_q = t( c( Mu1[q] - xi ,       Mu2[q] - delta))
    cov_q = make_cov(SigmaX[q],SigmaY[q], SigmaXY[q])
    
    covv_noise =  rho_xi*sigmaxi * sigmadelta
    cov_xi = make_cov( sigmaxi^2, sigmadelta^2 ,  covv_noise )
    
    
    SampleSigma = (cov_q + cov_xi )/2
    SampleSigma[is.nan(SampleSigma)] = 0
    
    test_stat = get_TwoSampleTestStat(dist_q, SampleSigma)
    
    test_stat
  }
  Test_XY_Noise <-function (q,   Mu1,Mu2,  SigmaX,SigmaY , SigmaXY,    xi, delta,sigmaxi, sigmadelta, rho_xi) {
    
    dist_q = t( c( Mu1[q] - xi ,       Mu2[q] - delta))
    cov_q = make_cov(SigmaX[q],SigmaY[q], SigmaXY[q])
    
    covv_noise =  rho_xi*sigmaxi * sigmadelta
    cov_xi = make_cov( sigmaxi^2, sigmadelta^2 ,  covv_noise )
    
    
    SampleSigma = (cov_xi  )
    SampleSigma[is.nan(SampleSigma)] = 0
    
    test_stat = get_mahala_ginv (dist_q, SampleSigma)
    
    test_stat
  }
}
