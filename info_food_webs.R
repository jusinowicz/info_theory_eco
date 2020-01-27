#=============================================================================
# R code to create to explore the Information Theoretic properties of 
# simple food webs. This creates a simple food web with an underlying dynamic
# model. 	  
# 1. Food-web includes resource, herbivore, and predator: 
#	 A. Resource is based on a consumer-resource model, with added predators 
#		1. Competition between consumers and resources emerges from consumption
#		2. Parameters at each level can be made a function of temperature. 
#	 B. Resources can be stochastic due to environmental fluctuations. 
#	 C. Relative non-linearity allows 2 consumers per Resource
# 2. Generate a bunch of random food webs 
# 3. Use information theory to track the resulting food-web structures. 
#=============================================================================
#=============================================================================
# load libraries
#=============================================================================
library(deSolve)
library(fields)
source("./food_web_functions.R")
source("./info_theory_functions.R")

#=============================================================================
# Outer loop. Set the number of trials and determine how to generate 
# combinations of species and parameters. 
#=============================================================================

#Length and time steps of each model run
tend = 200
delta1 = 0.01
tl=tend/delta1

#The maximum block depth for dynamic info metrics (larger is more accurate, but
#slower and could cause crashing if too large)
k= 10 

#Number of food webs to generate
nwebs = 1
#Output of each web
out1 = list(matrix(0,nwebs,1))
#Converting the web to Rutledge's compartment model and calculating the information
#theoretic quantities: Shannon Entropy, Mutual Information, Conditional Entropy
rweb1 = list(matrix(0,nwebs,1))
#Dynamic information metrics calculated from the (discretized) time series 
di_web = list(matrix(0,nwebs,1))

#Random resources:
c = 0.1
amp = 1
res_R = c(amp,c)

for (w in 1:nwebs){ 
	print(w)

	#Assume 3 trophic levels unless otherwise specified.
	nRsp = ceiling(runif(1)*30)
	nCsp = ceiling(runif(1)*20)
	nPsp = ceiling(runif(1)*10)
	nspp = nRsp+nCsp+nPsp

	#Randomly generate the species parameters for the model as well: 
	spp_prms = NULL
	#Resource: Nearly identical resource dynamics: 
	spp_prms$rR = matrix(rnorm(nRsp,300,10), nRsp, 1) #intrinsic growth
	spp_prms$Ki = matrix(rnorm(nRsp,500,10), nRsp, 1) #carrying capacity

	#Consumers: 
	spp_prms$rC = matrix(rnorm(nCsp,.5,0.2), nCsp, 1) #intrisic growth
	spp_prms$eFc = matrix(1,nCsp,nRsp) # just make the efficiency for everything 1 for now
	spp_prms$muC = matrix(rnorm(nCsp,0.6,0.1), nCsp, 1) #mortality rates
	#Consumption rates: 
	#Generate a hierarchy where each species predominantly feeds on particular resource. 
	dspp = abs((nCsp - nRsp))
	hier1= seq(1/nRsp, (1-1/nRsp), length=nRsp)
	spp_prms$cC = hier1 
	for( n in 1:nCsp) {
		spp_prms$cC = cbind(spp_prms$cC, shifter(hier1,n))
	}
	spp_prms$cC = matrix(spp_prms$cC[1:nRsp,1:nCsp ],nRsp,nCsp)

	#Predators: 
	spp_prms$rP = matrix(rnorm(nPsp,0.5,0.2), nPsp, 1) #intrisic growth
	spp_prms$eFp = matrix(1,nPsp,nCsp) # just make the efficiency for everything 1 for now
	spp_prms$muP = matrix(rnorm(nPsp,0.6,0.1), nPsp, 1) #mortality rates
	#Consumption rates: 
	#Generate a hierarchy where each species predominantly feeds on particular resource. 
	dspp = ((nPsp - nCsp))
	if(dspp<0){dspp = 0 }
	hier1= seq(1/nCsp, (1-1/nCsp), length = nCsp)
	spp_prms$cP = hier1
	for( n in 1:nPsp) {
		spp_prms$cP = cbind(spp_prms$cP, shifter(hier1,n))
	}
	spp_prms$cP = matrix(spp_prms$cP[1:nCsp,1:nPsp],nCsp,nPsp)


	#=============================================================================
	# Inner loop. Run the food web model, calculate information theoretic 
	# quantities. 
	#=============================================================================
	#=============================================================================
	# This function gives: 
	# out 		The time series for of population growth for each species in the web
	#			This can be set to just give the final 2 time steps of the web with
	#			"final = TRUE"
	# spp_prms	The parameters of all species in the food web
	#=============================================================================
	# tryCatch( {out1[w] = list(food_web_dynamics (spp_list = c(nRsp,nCsp,nPsp), spp_prms = spp_prms, 
	# 	tend, delta1, res_R = NULL,final = FALSE ))}, error = function(e){}) 
	#Random resource fluctuations:
	tryCatch( {out1[w] = list(food_web_dynamics (spp_list = c(nRsp,nCsp,nPsp), spp_prms = spp_prms, 
		tend, delta1, res_R = res_R,final = FALSE ))}, error = function(e){}) 

	
	#=============================================================================
	# Information theoretic assessment of the foodweb.
	#=============================================================================
	#=============================================================================
	# This section is as per Rutledge, Basore, and Mulholland 1976
	#=============================================================================
	## This code takes the ODEs and converts them to a biomass balance matrix and 
	## transition matrix. 
	## This version creates a compartment for each "event" where biomass is gained 
	## or loss. This includes birth, death, and "inefficiency" in the form of the 
	## way that biomass consumed translates to new population biomass. 
	#=============================================================================
	# This function gives:
	# Qi(t)		Biomass proportion flow through a node at time t
	# fij(t)	Probability of biomass flow between i and j at t
	# fijQi(t)  Total biomass flowing from i to j at t
	# sD 		Shannon entropy
	# mI_mean	Average mutual information
	# mI_per	Mutual information per interaction
	# ce 		Conditional entropy		
	#=============================================================================
	
	rweb1[w] = list(rutledge_web( spp_list=c(nRsp,nCsp,nPsp), pop_ts = out1[[w]]$out[,2:(nspp+1)],
		spp_prms = out1[[w]]$spp_prms) )

	#=============================================================================
	# Information processing networks
	#=============================================================================
	## This code takes the population time-series counts output by the ODEs and 
	## calculates Excess Entropy, Active Information Storage, and Transfer Entropy. 
	#=============================================================================
	# This function gives:
	# EE_mean		Average mutual information per species
	# AI_mean		Average active information per species
	# TE_mean		Average transfer entropy per species
	# 
	# EE_local		Local mutual information per species
	# AI_local		Local active information per species
	# TE_local		Local transfer entropy per species
	#=============================================================================
	nt1 = 1
	nt2 = tl
	di_web[w] = list(get_info_dynamics(pop_ts = floor(out1[[w]]$out[nt1:tl,2:(nspp+1)]), k=k ))

}

#=============================================================================
# Basic plots of species populations. 
#=============================================================================
w=1
out = out1[[w]]$out
nspp = out1[[w]]$spp_prms$nspp
nRsp = out1[[w]]$spp_prms$nRsp
nCsp = out1[[w]]$spp_prms$nCsp
nPsp = out1[[w]]$spp_prms$nPsp
tl = tend/delta1
plot(out[,"1"],t="l",col="red",ylim = c(0,max(out[tl,2:nspp],na.rm=T)))
for( n in 2:(nRsp) ) {
lines(out[,paste(n)],t="l",col="red")
}
for( n in ( (nRsp+1):(nRsp+nCsp) ) ) {
lines(out[,paste(n)],t="l",col="blue")
}
for( n in ((nRsp+nCsp):nspp ) ) {
lines(out[,paste(n)],t="l")
}

#=============================================================================
# Plot the dynamic information metrics with time 
#=============================================================================
#Local excess entropy
nt_use = dim(di_web[[w]]$ee_local)[1]
image.plot( 1:nt_use, 1:nspp, di_web[[w]]$ee_local, ylab="Species number", xlab="Time" )
abline(h =out1[[w]]$spp_prms$nRsp )
mtext("Resources", side=2, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(h =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp  )
mtext("Consumers", side=2, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Predators", side=2, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )

#Local active information storage
nt_use = dim(di_web[[w]]$ai_local)[1]
image.plot( 1:nt_use, 1:nspp, di_web[[w]]$ai_local, ylab="Species number", xlab="Time" )
abline(h =out1[[w]]$spp_prms$nRsp )
mtext("Resources", side=2, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(h =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp  )
mtext("Consumers", side=2, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Predators", side=2, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )

#Local transfer entropy
nt_use = dim(di_web[[w]]$te_local)[1]
image.plot( 1:nt_use, 1:nspp, di_web[[w]]$te_local, ylab="Species number", xlab="Time" )
abline(h =out1[[w]]$spp_prms$nRsp )
mtext("Resources", side=2, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(h =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp  )
mtext("Consumers", side=2, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Predators", side=2, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )

#Local separable information
nt_use = dim(di_web[[w]]$si_local)[1]
image.plot( 1:nt_use, 1:nspp, di_web[[w]]$si_local, ylab="Species number", xlab="Time" )
abline(h =out1[[w]]$spp_prms$nRsp )
mtext("Resources", side=2, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(h =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp  )
mtext("Consumers", side=2, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Predators", side=2, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )


###Or plot a subset of the data: 
nt1 = 5000
nt2 = tl-50
image.plot( nt1:nt2, 1:nspp, di_web[[w]]$ee_local[nt1:nt2,], ylab="Species number", xlab="Time" )

#Local active information storage
image.plot( nt1:nt2, 1:nspp, di_web[[w]]$ai_local[nt1:nt2,], ylab="Species number", xlab="Time" )

#Local transfer entropy
image.plot( nt1:nt2, 1:nspp, di_web[[w]]$te_local[nt1:nt2,], ylab="Species number", xlab="Time" )

#Local separable information
image.plot( nt1:nt2, 1:nspp, di_web[[w]]$si_local[nt1:nt2,], ylab="Species number", xlab="Time" )
abline(h =out1[[w]]$spp_prms$nRsp )
mtext("Resources", side=2, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(h =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp  )
mtext("Consumers", side=2, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Predators", side=2, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )



out1[[w]]$out[10000,]>1e-5

#Generate quantities for the maximum entropy distribution, i.e. uniform: 
pop_me = runif(nspp)
me_freq = pop_me/matrix(sum(pop_me),length(pop_me),1)

