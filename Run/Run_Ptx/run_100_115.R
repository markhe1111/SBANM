setwd('/Users/markhe 1/Dropbox/DifferentialBlockModel/Github/')

rm(list = ls()) 
load('Data/Politics/REproc_100_115.Rdata')

# load packages
library(matlib)
library(DescTools)   
library(Brobdingnag)
library(igraph) 
library(MASS)
library(fpc)
library(rlist) 
library(randnet)

Rcpp::sourceCpp('cpp/Bivariate/E_Biv.cpp')
Rcpp::sourceCpp('cpp/Bivariate/MStepBiv.cpp')
Rcpp::sourceCpp('cpp/Bivariate/biv_ELBO.cpp')
#Rcpp::sourceCpp('cpp/ClosedP_q/M_Step_b.cpp')

source('SBANM/Utilities/general_Utils.R')
source('SBANM/Utilities/Biv_Utils.R')
source('SBANM/Bivariate.R')

A1=  FisherZ( adjust(  AA2,  .0001)  )
A2=  FisherZ( adjust(  AA3, .0001)  )
 
diag(A1) = diag(A2) = 0
   
# running example 

P3Gps = bivSBANM (Q=3, A1,  A2,   max_iter = 200, StochSize.a = 100 , TOL =1e-12 )
P4Gps = bivSBANM (Q=4, A1,  A2,   max_iter = 200, StochSize.a = 100 , TOL =1e-12 )
