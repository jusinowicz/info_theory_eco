#=============================================================================
# R code to measure the fitness value of information for a deterministic 
# model of competitive community dynamics. This is to first attempt, to ask 
# the  question, is it even meaningful to think in these terms for a 
# deterministic model? 
# 
# Model from Gil et al. 2019 
#=============================================================================
#=============================================================================
#Load libraries
#=============================================================================
library(tidyverse)
library(RandomFields)
library(vegan)
library(deSolve)
source("./env_functions.R")

#=============================================================================
# Define the population dynamics through the following functions
#=============================================================================
#No information: 
rc_noi = function(times,sp,parms){
	with( as.list(c(parms, sp )),
		{	
    P = matrix(sp[1:nPsp],nPsp, 1)

	###Resource (prey) dynamics: Logistic growth, reduced by consumption
	dP = P
	for( i in 1:nspp){
	dP[i] = P[i]*( (rp[i]) * (1 - alphas[i,] %*%P/Kp[i]) - B[i])
	}
	# ###Consumer dynamics: LV consumption
	# dH = H 
	# dH = H*( ah*P - H*mh )	

	list(c(dP)) 
	})	
}
  
#With social information: 
rc_i = function(times,sp,parms){
	with( as.list(c(parms, sp )),
		{	
    P = matrix(sp[1:nPsp],nPsp, 1)

	###Resource (prey) dynamics: Logistic growth, reduced by consumption
	dP = P
	for( i in 1:nspp){
	dP[i] = P[i*( (rp[i]) * (1 - alphas[i,] %*%P/Kp[i]) - ps[i]*(exp( -(b[i,]%*%P) ) )-pm )
	}

	# ###Consumer dynamics: LV consumption
	# dH = H 
	# dH = H*( ah*P - H*mh )	

	list(c(dP)) 
	})	
}
  
#==============================================================================
# Set parameter values 
#==============================================================================
###Time to simulate over: 
tend = 1000 #Length of numerical simulation
times  = seq(from = 0, to = tend, by = 0.01)
tl = length(times)

####No information 
spp_prms_noi = NULL
spp_prms_noi$nPsp = 1 # number of species 
spp_prms_noi$rp = 8 #intrinsic growth
spp_prms_noi$Kp = 100 #carrying capacity
spp_prms_noi$B = 5 #predation rate

####With social information
spp_prms_i = NULL
spp_prms_i$nPsp = 1 # number of species 
spp_prms_i$rp = 8 #intrinsic growth
spp_prms_i$Kp = 100 #carrying capacity
spp_prms_i$ps = 3 #effect of social information
spp_prms_i$pm = 2 #min predation rate
spp_prms_i$b = 0.1 #social information rate


parms_noi = list(nPsp=spp_prms_noi$nPsp, rp = spp_prms_noi$rp,
	Kp =spp_prms_noi$Kp, B = spp_prms_noi$B 
 )

parms_i = list(nPsp=spp_prms_i$nPsp, rp = spp_prms_i$rp,
	Kp =spp_prms_i$Kp, ps = spp_prms_i$ps, pm = spp_prms_i$pm, 
	b = spp_prms_i$b 
 )
#=============================================================================
# Run the model with initial conditions
#=============================================================================
minit = c( matrix(25,spp_prms$nPsp,1) )
####No information 
rc_noi_out = ode(y=minit, times=times, func=rc_noi, parms=parms_noi, atol = 1e-9)
rc_noi_out = as.data.frame(rc_noi_out)

####With social information
rc_i_out = ode(y=minit, times=times, func=rc_i, parms=parms_i, atol = 1e-9)
rc_i_out = as.data.frame(rc_i_out)

#=============================================================================
# Plot
#=============================================================================
#Change to long format: 
rc_noi_long =rc_noi_out %>% gather( species, N, 2:(spp_prms$nPsp+1) )
rc_i_long =rc_i_out %>% gather( species, N, 2:(spp_prms$nPsp+1) )

#Merge data frames: 
rc_i_long$species[rc_i_long$species=="1"] = "2"
both_long = rbind(rc_noi_long ,rc_i_long ) #inner join

#Each species' time trajectory
p0 = ggplot()+ geom_line( data = both_long, aes ( x = time, y = N, color = species)  )+ 
ylab("Population")+
theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position = "none") 



####
#Each species' time trajectory
p1=ggplot()+ geom_line( data = rc_noi_long, aes ( x = time, y = N, color = species)  )+ 
ylab("Population")+
theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position = "none") 

#########
rc_i_long =rc_i_out %>% gather( species, N, 2:(spp_prms$nPsp+1) )

#Each species' time trajectory
p2=ggplot()+ geom_line( data = rc_i_long, aes ( x = time, y = N, color = species)  )+ 
ylab("Population")+
theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position = "none") 


#p4 = grid.arrange(p1,p2 nrow = 2 )

