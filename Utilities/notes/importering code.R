
TX = makeTX(tau, A1)
TY = makeTX(tau, A2)
TZ = makeTX(tau, A3)
T0 = makeT0(tau )

TX_Q = makeTX_PQ (tau, A1, P.00 )
TY_Q = makeTX_PQ (tau, A2, P.00)
TZ_Q = makeTX_PQ (tau, A3, P.00)

T0_Q = makeT0_Q(tau,  P.00) 

xi1 =    nan_omit(TX / T0)   
xi2 =    nan_omit(TY / T0)    
xi3 =     nan_omit(TZ / T0)    
xis = c(xi1,xi2,xi3)

xiB1 =     nan_omit( TX_Q /  T0_Q )
xiB2 =     nan_omit( TY_Q /  T0_Q )
xiB3 =    nan_omit( TZ_Q /  T0_Q )
xiBs = c(xiB1,xiB2,xiB3)

Mu1 = MStep_Mu_Q(A1,tau) * P.11 + (1-P.11)* xiB1
Mu2 = MStep_Mu_Q(A2,tau  )* P.11+ (1-P.11)* xiB2
Mu3 = MStep_Mu_Q(A3,tau  )* P.11+ (1-P.11)* xiB3

TX2 =  makeTX2(A1,tau,xi1)
TY2 =  makeTX2(A2,tau,xi2)
TZ2 =  makeTX2(A3,tau,xi3)

# PQ is in alt
TX2Q = makeTX2_PQ  (tau,A1, P00 = P.00   ,xiB1)
TY2Q = makeTX2_PQ  (tau,A2, P00 = P.00   ,xiB2)
TZ2Q = makeTX2_PQ  (tau,A3, P00 = P.00   ,xiB3)

Sigmaxi1  =   nan_omit(TX2 /T0) # * Psi
Sigmaxi2  =   nan_omit(TY2 /T0)   
Sigmaxi3  =   nan_omit(TZ2 /T0 )  

SigmaxiB1  =     nan_omit(TX2Q /  T0_Q ) #*(1-Psi)
SigmaxiB2  =     nan_omit(TY2Q /  T0_Q )
SigmaxiB3  =     nan_omit(TZ2Q /  T0_Q )

sigmaXis = c(sqrt(Sigmaxi1),sqrt(Sigmaxi2),sqrt(Sigmaxi3))
sigmaXiBs = c(sqrt(SigmaxiB1),sqrt(SigmaxiB2),sqrt(SigmaxiB3))

SigmaX = MStep_SigmaX ( A1,tau,Mu1, xi1)* P.11 + SigmaxiB1*P.00
SigmaY = MStep_SigmaX ( A2,tau,Mu2, xi2)* P.11 + SigmaxiB2*P.00
SigmaZ = MStep_SigmaX ( A3,tau,Mu3, xi3)* P.11 + SigmaxiB3*P.00

SigmaXY = MStep_SigmaXY (A1,A2,tau,Mu1,Mu2 ,xiB1,xiB2 )*P.11
SigmaYZ = MStep_SigmaXY (A2,A3,tau,Mu2,Mu3 ,xiB2,xiB3 )*P.11
SigmaXZ = MStep_SigmaXY (A1,A3,tau,Mu1,Mu3 ,xiB1,xiB3 )*P.11

rhos12 = SigmaXY  / sqrt(SigmaX * SigmaY)
rhos23 = SigmaYZ  / sqrt(SigmaY * SigmaZ)
rhos13 = SigmaXZ  / sqrt(SigmaX * SigmaZ)