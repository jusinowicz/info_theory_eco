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
library(gridExtra)
library(viridis)
source("./env_functions.R")

#=============================================================================
# Define the population dynamics through the following functions
#=============================================================================
#No information: 
cc_noi = function(times,sp,parms){
	with( as.list(c(parms, sp )),
		{	
    P = matrix(sp[1:nPsp],nPsp, 1)

	###Resource (prey) dynamics: Logistic growth, reduced by consumption
	dP = P
	for( i in 1:nspp){
	dP[i] = P[i]*( rp[i]* ( (Kp[i] - (alphas[i,]%*%P ) )/Kp[i] ) - B[i] )
	}
	# ###Consumer dynamics: LV consumption
	# dH = H 
	# dH = H*( ah*P - H*mh )	

	list(c(dP)) 
	})	
}
  
#With social information: 
cc_i = function(times,sp,parms){
	with( as.list(c(parms, sp )),
		{	
    P = matrix(sp[1:nPsp],nPsp, 1)

	###Resource (prey) dynamics: Logistic growth, reduced by consumption
	dP = P
	for( i in 1:nspp){
	dP[i] = P[i]*( rp[i] *( (Kp[i]- (alphas[i,] %*%P) )/Kp[i] )- (ps[i]*exp( -(b[i,]%*%P) )+pm[i]) )
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
spp_prms_noi$nPsp = 2 # number of species 
spp_prms_noi$rp = c(5,5) #intrinsic growth
spp_prms_noi$Kp = c(100,100) #carrying capacity
spp_prms_noi$B = c(4,4) #predation rate
spp_prms_noi$alphas = 
		matrix( c(1,.8,.8,1), #c(1,1.5,1.5,1)
		spp_prms_noi$nPsp,spp_prms_noi$nPsp) #competition


####With social information
spp_prms_i = NULL
spp_prms_i$nPsp = 2 # number of species 
spp_prms_i$rp = c(5,5) #intrinsic growth
spp_prms_i$Kp = c(100,100) #carrying capacity
spp_prms_i$ps = c(2,2) #effect of social information
spp_prms_i$pm = c(2,2) #min predation rate
spp_prms_i$alphas = 
		matrix( c(1,.8,.8,1), #c(1,1.5,1.5,1)
		spp_prms_i$nPsp,spp_prms_i$nPsp) #competition
spp_prms_i$b = 
		matrix( c(0,0.1,0.1,0),
		spp_prms_i$nPsp,spp_prms_i$nPsp) #information rate


parms_noi = list(nPsp=spp_prms_noi$nPsp, rp = spp_prms_noi$rp,
	Kp =spp_prms_noi$Kp, B = spp_prms_noi$B, alphas = spp_prms_noi$alphas
 )

parms_i = list(nPsp=spp_prms_i$nPsp, rp = spp_prms_i$rp,
	Kp =spp_prms_i$Kp, ps = spp_prms_i$ps, pm = spp_prms_i$pm, 
	b = spp_prms_i$b, alphas = spp_prms_i$alphas
 )

#=============================================================================
# Run the model with initial conditions
#=============================================================================
minit = c( matrix( c(.01,20),spp_prms_i$nPsp,1) )
####No information 
cc_noi_out = ode(y=minit, times=times, func=cc_noi, parms=parms_noi, atol = 1e-9)
cc_noi_out = as.data.frame(cc_noi_out)

####With social information
cc_i_out = ode(y=minit, times=times, func=cc_i, parms=parms_i, atol = 1e-9)
cc_i_out = as.data.frame(cc_i_out)

#=============================================================================
# Plot
#=============================================================================
#Change to long format: 
cc_noi_long =cc_noi_out %>% gather( species, N, 2:(spp_prms_noi$nPsp+1) )
cc_i_long =cc_i_out %>% gather( species, N, 2:(spp_prms_i$nPsp+1) )

#Merge data frames: 
cc_i_long$species[cc_i_long$species=="1"] = "3"
cc_i_long$species[cc_i_long$species=="2"] = "4"
both_long = rbind(cc_noi_long ,cc_i_long ) #inner join

#Each species' time trajectory
p0 = ggplot()+ geom_line( data = both_long, aes ( x = time, y = N, color = species)  )+ 
ylab("Population")+
theme(axis.text.x=element_blank(), axis.title.x=element_blank()) #, legend.position = "none") 
p0

#=============================================================================
# Calculate the fitness value of information
#	This is done by subtracting the instantaneous boundary growth rate 
#	(invasion growth rate) of the no-info from that of the info model. 
#=============================================================================
nPsp = spp_prms_noi$nPsp
####No information
Pi = matrix(c(0,20,20,0),nPsp,nPsp) #Resident equilibrium population
igr_noi = spp_prms_noi$rp* ( (spp_prms_noi$Kp - 
			rowSums(spp_prms_noi$alphas*Pi ) )/spp_prms_noi$Kp) - 
			spp_prms_noi$B

####With social information
Pi2 = matrix( c(0,20,20,0 ),nPsp,nPsp)#Resident equilibrium population
igr_i =  spp_prms_i$rp*( (spp_prms_i$Kp- 
			rowSums(spp_prms_i$alphas*Pi2 ) )/spp_prms_i$Kp ) - 
			(spp_prms_i$ps*exp( -(rowSums(spp_prms_i$b*Pi2)) )+
			 spp_prms_i$pm)

####Fitness value of information 
fvoi = igr_i-igr_noi

#Which is analytically equivalent to: 
social_info = (spp_prms_i$ps*exp( -(rowSums(spp_prms_i$b*Pi2)) )+
			 spp_prms_i$pm)

fvoi2 = spp_prms_noi$B - social_info

#=============================================================================
# Paper plots: See figures.R
#=============================================================================

