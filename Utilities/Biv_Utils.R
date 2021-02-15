

make_variances_Biv  <-function(varX_qq, varY_qq,rho_qq){
  Q = length(varX_qq)
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
get_EN_each_q <-function(q, A1 ,   A2 , tau ,  Mu1,  Mu2, Xis,
                         SigmaX, SigmaY, SigmaXY, sigmaXis){
  
  CovMat_q = makeCovSig( q-1, SigmaX,SigmaY , SigmaXY)  
  covNoi = diag( sigmaXis^2) 
  Mus_q = cbind(Mu1,Mu2)[q,]
  Freq_q = E_Nq_Simp_Biv   ( A1,A2, tau[,(q)],
                             Mus_q, CovMat_q,  ginv(CovMat_q),  
                             Xis,   covNoi,    ginv(covNoi))
  Freq_q
}

get_bivELBO <- function(A1,A2, try){
  elbo_values = bivELBO ( A1,A2,  try$tau  ,     
                          try$Mu1,  try$Mu2,   try$SigmaX, try$SigmaY, 
                          try$SigmaXY,  try$Xis ,   try$sigmaXis  ,   try$P.11)
  sum( elbo_values  [ upper.tri(elbo_values) ])
  
}

make_cov2 <-function(SigmaX_i, SigmaY_i, SigmaXY_i ){
  covar = matrix(ncol=2, nrow=2,0)
  covar[1,1] = SigmaX_i
  covar[1,2] = covar[2,1] = SigmaXY_i
  covar[2,2] = SigmaY_i
  covar
}

get_EN_each_q <-function(q, A1 ,   A2 , tau , 
                         Mu1,  Mu2, Xis,
                         SigmaX, SigmaY, SigmaXY, sigmaXis){
  
  CovMat_q = makeCovSig( q-1, SigmaX,SigmaY , SigmaXY)  
  covNoi = diag( sigmaXis^2) 
  Mus_q = cbind(Mu1,Mu2)[q,]
  Freq_q = E_Nq_Simp_Biv   ( A1,A2, tau[,(q)],
                             Mus_q, CovMat_q,  ginv(CovMat_q),  
                             Xis,   covNoi,    ginv(covNoi))
  Freq_q
}



get_bivELBO <- function(A1,A2, try){
  elbo_values = bivELBO ( A1,A2,  try$tau  ,     
                          try$Mu1,  try$Mu2,   try$SigmaX, try$SigmaY, 
                          try$SigmaXY,  try$Xis ,   try$sigmaXis  ,   try$P.11)
  sum( elbo_values  [ upper.tri(elbo_values) ])
  
}
get_ICL2 = function(A1,A2 ,try , Q,N ){
  
  ELBO =     get_bivELBO( A1,A2,  try ) 
  pen = Q*log(N*(N-1)*2/2) + Q*(Q-1)/2 * 2 * log(N*(N-1)/2)
  ICL = ELBO -.5*Q*(Q-1) *log(N) - pen
  ICL
}
