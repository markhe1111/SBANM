make_var_Triv  <-function(Q, varX_qq, varY_qq, varZ_qq, rhos){
  
  colnames(rhos) = c('cor_anx_beh', 'cor_beh_mood', 'cor_anx_mood')
  #rownames(rhos ) = c('B1','B2','B3')
  Rhos = data.frame(rhos)
  
  variances = list()
  
  for(q in 1:Q){
    varmat = matrix(ncol=3, nrow=3, 0)
    colnames(varmat) =  c('anx', 'beh', 'mood')
    rownames(varmat) =  c('anx', 'beh', 'mood')
    
    varmat[1,1] = varX_qq[q]
    varmat[2,2] = varY_qq[q]
    varmat[3,3] = varZ_qq[q]
    
    varmat[1,2] = sqrt(varX_qq[q] * varY_qq[q]) *  Rhos$cor_anx_beh[q]
    varmat[2,1] = sqrt(varX_qq[q] * varY_qq[q]) * Rhos$cor_anx_beh[q]
    
    varmat[1,3] = sqrt(varX_qq[q] * varZ_qq[q]) * Rhos$cor_anx_mood[q]
    varmat[3,1] = sqrt(varX_qq[q] * varZ_qq[q]) * Rhos$cor_anx_mood[q]
    
    varmat[2,3] = sqrt(varY_qq[q] * varZ_qq[q]) * Rhos$cor_beh_mood[q]
    varmat[3,2] = sqrt(varY_qq[q] * varZ_qq[q]) * Rhos$cor_beh_mood[q]
    
    variances[[q]] = varmat
  }
  return(variances)
}

make_variances_Trivariate  <-function(Q,varX_qq, varY_qq, varZ_qq,rho_qq){
  variances = list()
  for(q in 1:Q){
    varmat = matrix(ncol=3, nrow=3, 0)
    
    varmat[1,1] = varX_qq[q]
    varmat[2,2] = varY_qq[q]
    varmat[3,3] = varZ_qq[q]
    
    varmat[1,2] = sqrt(varX_qq[q] * varY_qq[q]) * rho_qq[q]
    varmat[2,1] = sqrt(varX_qq[q] * varY_qq[q]) * rho_qq[q]
    
    varmat[1,3] = sqrt(varX_qq[q] * varZ_qq[q]) * rho_qq[q]
    varmat[3,1] = sqrt(varX_qq[q] * varZ_qq[q]) * rho_qq[q]
    
    varmat[2,3] = sqrt(varY_qq[q] * varZ_qq[q]) * rho_qq[q]
    varmat[3,2] = sqrt(varY_qq[q] * varZ_qq[q]) * rho_qq[q]
    
    variances[[q]] = varmat
  }
  return(variances)
}
  
 
make_SigmaCov <-function( q, SigmaX,SigmaY , SigmaZ,  SigmaXY, SigmaYZ ,  SigmaXZ){
  
  CovM = matrix(ncol= 3, nrow = 3)
  diag(CovM ) = c(SigmaX[q] ,SigmaY[q] , SigmaZ[q])
  CovM[1,2] =    CovM[2,1] =   SigmaXY[q]
  CovM[2,3] =    CovM[3,2] =   SigmaYZ[q]
  CovM[1,3] =    CovM[3,1] =   SigmaXZ[q]
  CovM
  
}
 
make_cor_mat  <-function(q , rhos12 , rhos23, rhos13){
  
  CorMat = diag(1,3)
  CorMat [1,2] =    CorMat [2,1] = rhos12 [q]
  CorMat [2,3] =    CorMat [3,2] = rhos23 [q]
  CorMat [1,3] =    CorMat [3,1] = rhos13 [q]
  CorMat
  
}



get_EN_each_q <-function(q, A1 ,   A2,  A3 , tau ,  
                         Mu1,  Mu2,  Mu3,Xis, 
                         SigmaX, SigmaY, SigmaZ,
                         SigmaXY, SigmaYZ,  SigmaXZ, sigmaXis){
  
  CovMat_q = makeCovSig( q-1, SigmaX,SigmaY , SigmaZ, SigmaXY, SigmaYZ,SigmaXZ)  
  covNoi = diag( sigmaXis^2) 
  Mus_q = cbind(Mu1,Mu2,Mu3)[q,]
  Freq_q = E_Nq_Simp_q   ( A1,A2, A3, tau[,(q)],
                           Mus_q, CovMat_q,  ginv(CovMat_q),  
                           Xis,   covNoi,    ginv(covNoi))
  Freq_q
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
