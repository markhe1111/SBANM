 
setwd('~/Dropbox/DifferentialBlockModel/multiSBM/data/Politics/other_files/')

load('ICL_find_min5.Rdata')

source('~/Dropbox/DifferentialBlockModel/multiSBM/util/Tri_Utils.R')

districts = colnames(A1)

setwd('/Users/markhe 1/Dropbox/DifferentialBlockModel/multiSBM/data/Politics/other_files/')
kos = read.csv('Daily Kos Elections urban-suburban-rural population by congressional district - Total pop..csv')

memberships = apply( try_from2_to7[[5]]$tau, 1, which.max)

MM1 = M2_ab
MM2 = M3_ab
Q = 5
RUN =  try_from2_to7[[5]]

P. = RUN$P.11
tautable(RUN$tau)

memQ = sort(unique(memberships))
gps = lapply(1:Q, function(q) which(memberships==q))
blocs = lapply(gps, function(g) districts[g])


B1_all  =    lapply(1:Q, function(q)   MM1  [ as.character(MM1$distcode) %in% blocs[[q]],] )
B2_all  =    lapply(1:Q, function(q)   MM2  [ as.character(MM2$distcode) %in% blocs[[q]],] )


get_avg_repub = function(Block){
  mean ( Block$party_code==200 )
}
 
Avg_repub = cbind( sapply(B1_all, function(x)  get_avg_repub(x) ),
                   sapply(B2_all, function(x)  get_avg_repub(x) ) )
 
#  table(memberships)
rhos =RUN$rhos 
RUN$P.11
Mu1 = RUN$Mu1
Mu2 = RUN$Mu2

gopmu_table = cbind(   table(memberships), Mu1, Mu2 , round(rhos ,2),   round( Avg_repub,2))


colnames(gopmu_table) = c('n','mu_1','mu_2', 'rho','GOP100','GOP115')

xtable::xtable(gopmu_table)

B1_all

B1_all[[1]]
B2_all[[4]]


dont_worry_abt_kos_stuff <-function(){

  Avg_age = cbind( sapply(B1_all, function(x)    get_avg_age_100(x) ),
                   sapply(B2_all, function(x)    get_avg_age_115(x) ) )
  
  
  get_avg_age_100 = function(Block){   mean( 1987- Block$born )  }
  get_avg_age_115 = function(Block){    mean ( 2017- Block$born )  }
  
# kos stuff

kos_state = substr(as.character( kos$District),0,2)
kos_dist = substr(as.character( kos$District),4,5)
as.numeric(kos_dist)
kos_state_dist = paste(kos_state,   as.numeric(kos_dist), sep = ' ')
kos$dist_code = kos_state_dist
B_kos  =    lapply(1:Q, function(q)   kos  [  kos$dist_code %in% blocs[[q]] ,  ] )

get_avg_urban = function(kosBlock){
  urban = mean ( kosBlock$Urban.. )
  suburban = mean ( kosBlock$Suburban.. )
  rural = mean ( kosBlock$Rural.. )
  c(urban, suburban, rural)
}

get_urban = function(kosBlock){
  mean ( kosBlock$Most.Prevalent=='Urban')
}
get_rural = function(kosBlock){
  mean ( kosBlock$Most.Prevalent=='Rural')
}
get_subUrban = function(kosBlock){
  mean ( kosBlock$Most.Prevalent=='Suburban')
}


Avg_urbtype =  sapply(B_kos, function(x)  get_avg_urban(x) )
Avg_urban =  sapply(B_kos, function(x)  get_urban(x) )
Avg_rural =  sapply(B_kos, function(x)  get_rural(x) )
Avg_suburban =  sapply(B_kos, function(x)  get_subUrban(x) )

 gopurb_table = cbind( gop_table , round( Avg_urban,2) ,   round( Avg_rural,2) ,   round( Avg_suburban,2))


colnames(gopurb_table) = c('n','rho','GOP100','GOP115','urb','rural')

xtable::xtable(gopurb_table)





param_table  = data.frame(noiseProb =P.[,2], 
                          Mu1 ,sigmax=sqrt(SigmaX), 
                          Mu2, sigmay=sqrt(SigmaY),
                          Mu3, sigmaz=sqrt(SigmaZ),
                          rhos12 ,rhos23,rhos13,
                          n=tautable_fill(tau, Q)) 


  #born
  # nominate_dim1 ,nominate_dim2, nominate_log_likelihood,
  #  nominate_geo_mean_probability, nominate_number_of_votes,
  #  nominate_number_of_errors,
  #  nokken_poole_dim1, nokken_poole_dim2)
  #nominate_dim1 ,nominate_dim2, nominate_log_likelihood,
  #nominate_geo_mean_probability, nominate_number_of_votes,
  #nominate_number_of_errors,
  #nokken_poole_dim1, nokken_poole_dim2)
  
  #View( cbind(colnames(A1) , colnames(A2)))
}
   
