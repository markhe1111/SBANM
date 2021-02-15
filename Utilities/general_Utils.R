adjust = function(X, epsilon){
  X[X==1] = 1 - epsilon
  X[X== -1] = -1 + epsilon
  X
} 
rowMax = function(M){
  apply(M, 1, max)
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
