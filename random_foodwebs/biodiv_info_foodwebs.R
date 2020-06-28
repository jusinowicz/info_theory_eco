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
# 4. This file has a lot of code for visualizing output of both the foodweb 
#	 its information theoretic properties after the main loop. 
#=============================================================================
#=============================================================================
# load libraries
#=============================================================================
library(deSolve)
library(fields)
library(tidyverse)
library(lubridate)
library(mgcv)
source("../info_theory_functions/food_web_functions.R")
source("../info_theory_functions/info_theory_functions.R")


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
k= 5 

#Number of food webs to generate
nwebs = 50

#Output of each web
out1 = vector("list",nwebs)
#Converting the web to Rutledge's compartment model and calculating the information
#theoretic quantities: Shannon Entropy, Mutual Information, Conditional Entropy
rweb1 = vector("list",nwebs)
#Dynamic information metrics calculated from the (discretized) time series 
di_web = vector("list",nwebs)
#Track the average transfer entropy and separable information between each pair of 
#species as a way to build a network of information flow through the network. 
te_web = vector("list",nwebs)
si_web = vector("list",nwebs)
aiE_web = vector("list",nwebs)
MMI_web = vector("list",nwebs)

#Random resources:
c = 0.1
amp = 1
res_R = c(amp,c)

for (w in 1:nwebs){ 
	print(w)

	#Assume 3 trophic levels unless otherwise specified.
	nRsp = ceiling(runif(1)*5)
	nCsp = ceiling(runif(1)*4)
	nPsp = ceiling(runif(1)*3)
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
		tend, delta1, res_R = res_R) )

		# print( paste( "nRsp", sum(out1[[w]]$out[tl,1:nRsp]>1) ) )
		# print( paste( "nCsp", sum(out1[[w]]$out[tl,(nRsp+1):nCsp]>1) ) )
		# print( paste( "nPsp", sum(out1[[w]]$out[tl,(nCsp+1):nPsp]>1) ) )		

		# plot(out1[[w]]$out[,1], t="l", ylim = c(0, max(out1[[w]]$out[tl,],na.rm=T) ) )
		# for(n in 2:nRsp){ lines(out1[[w]]$out[,n], col ="red") }
		# for(n in (nRsp+1):(nCsp) ){ lines(out1[[w]]$out[,n], col ="blue") }
		# for(n in (nCsp+1):(nPsp) ){ lines(out1[[w]]$out[,n]) }

	
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
	## Each quantity is calculated at both the average and local level.  
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
	# di_web[w] = list(get_info_dynamics(pop_ts = floor(out1[[w]]$out[nt1:tl,2:(nspp+1)]), 
	# 	k=k,with_blocks=FALSE))

	# ## This code takes the population time-series counts output by the ODEs and 
	# ## calculates the average Transfer Entropy from each species to every other 
	# ## species. The goal is to get an overview of the major information pathways 
	# ## in the web.   
	# #=============================================================================
	# # This function gives:
	# # te_web		Average transfer entropy per species as a pairwise matrix
	# #=============================================================================
	# te_web[w] = list( get_te_web( pop_ts = floor(out1[[w]]$out[nt1:tl,2:(nspp+1)]), 
	# 	k=k) )

	# ## This code takes the population time-series counts output by the ODEs and 
	# ## calculates the average Separable Information from each species to every other 
	# ## species. The goal is to get an overview of the major information pathways 
	# ## in the web.   
	# #=============================================================================
	# # This function gives:
	# # si_web		Average separable information per species as a pairwise matrix
	# #=============================================================================
	# si_web[w] = list( get_si_web( pop_ts = floor(out1[[w]]$out[nt1:tl,2:(nspp+1)]), 
	# 	k=k) )

	#=============================================================================
	# This function gives:
	# aiE_web    The AI of the entire ensemble, treated as a single time series. 
	#=============================================================================
	aiE_web[w] = list( get_ais (  series1 = floor(out1[[w]]$out[nt1:tl,2:(nspp+1)]), 
		k=k, ensemble = TRUE)    )

	#=============================================================================
	# This function gives:
	# MMI_web    The MMI of the entire ensemble, treated as a single time series. 
	#=============================================================================
	MMI_web[w] = list( get_ais (  series1 = floor(out1[[w]]$out[nt1:tl,2:(nspp+1)]), 
		k=k, ensemble = TRUE)    )

	}, error = function(e){}) 

}

#save(file = "rand_fwebmod6F.var", out1,  di_web,te_web,si_web)
save(file = "rand_fwebmod7D.var", out1, rweb1,aiE_web,MMI_web)


#=============================================================================
# Load saved foodwebs and look at relationships between function and various 
# measures of complexity. 
#=============================================================================
#This requires user input!
variable.list=list("out1", "di_web", "te_web","si_web")


file.name.list=c(
  "rand_fwebmod7A.var", 
  "rand_fwebmod7B.var",  
  "rand_fwebmod7C.var",
  "rand_fwebmod7D.var",
  "rand_fwebmod7E.var",
  "rand_fwebmod7F.var"
  )


#Combine the variables from each scenario file into one variable
var.length=length(variable.list)
nscen = length(file.name.list)
out1_all=NULL
rweb1_all = NULL
aiE_web_all = NULL
MMI_web_all = NULL

for (g in 1:nscen){
	load(file.name.list[[g]])
	nwebs = length( aiE_web[(!sapply(aiE_web,is.null) ) ])
	rweb1_all=c(rweb1_all, rweb1[(!sapply(rweb1,is.null) ) ])
	aiE_web_all=c(aiE_web_all, aiE_web[(!sapply(aiE_web,is.null) ) ])
	MMI_web_all=c(MMI_web_all, MMI_web[(!sapply(MMI_web,is.null) ) ])
	out1_all=c(out1_all, out1[1:nwebs])

}

#Take variables out of the lists to plot: 
ncells=length(aiE_web_all)
rDIT = data.frame(matrix(0, nrow=ncells, ncol =8) ) 
ncnames = c("fwno","Biomass", "var_Biomass", "nspp", "shannon", "rMI", "MI", "AI" )
colnames(rDIT) = ncnames

for (n in 1:ncells){
	rDIT$fwno[n] = n 
	rDIT$nspp[n] = out1_all[[n]]$spp_prms$nspp #Number of species

	tlast = dim(out1_all[[n]]$out)[1] - 1 #Length of time series
	rDIT$Biomass[n] = sum(out1_all[[n]]$out[tlast, 2:(rDIT$nspp[n]+1) ]) #Biomass at last time

	tbck = tlast*3/4 #Use a subset that excludes transient stage for variance
	rDIT$var_Biomass[n] = var( rowSums( out1_all[[n]]$out[ (tlast-tbck):tlast, 2:(rDIT$nspp[n]+1) ]) )

	#Shannon Diversity
	pi = out1_all[[n]]$out[tlast, 2:(rDIT$nspp[n]+1) ] / rDIT$Biomass[n]
	pi[pi <= 0 ] = NA
	rDIT$shannon[n] = - sum( pi*log(pi),na.rm=T )

	#Rutledge MI 
	rDIT$rMI[n] = rweb1_all[[n]]$mI_mean2

	#Multiple Mutual Information
	rDIT$MI[n] = MMI_web_all[[n]]$mean

	#Ensemble active information
	rDIT$AI[n] = aiE_web_all[[n]]$mean

}

lm_nspp = lm(rDIT$Biomass~rDIT$nspp)
lm_H=lm(rDIT$Biomass~rDIT$shannon)
lm_rMI=lm(rDIT$Biomass~rDIT$rMI)
lm_MI=lm(rDIT$Biomass~rDIT$MI)
lm_AI=lm(rDIT$Biomass~rDIT$AI)

#Plots 
ggplot (rDIT, aes(x = nspp, y = Biomass,color = "1" ) ) + geom_point () + 
	geom_point (aes(x = shannon, y =Biomass,color = "2")) +
	geom_point( aes (x = rMI, y=Biomass,color = "3" ) ) +
	geom_point( aes (x = MI, y=Biomass,color = "4" ) ) +
	geom_point( aes (x = AI, y=Biomass,color = "5" ) ) +
	#scale_color_discrete(name ="", labels = c("# Species", "SDI", "rMI","MI","AI" ) )


#Files to load
# file.name.list=c(
#   "rand_fwebmod6A.var", 
#   "rand_fwebmod6B.var",  
#   "rand_fwebmod6C.var",
#   "rand_fwebmod6D.var",
#   "rand_fwebmod6E.var",
#   "rand_fwebmod6F.var"
#   )



# #Combine the variables from each scenario file into one variable
# var.length=length(variable.list)
# nscen = length(file.name.list)
# out1_all=NULL
# di_web_all = NULL
# te_web_all = NULL
# si_web_all = NULL

# for (g in 1:nscen){
# 	load(file.name.list[[g]])
# 	nwebs = length( di_web[(!sapply(di_web,is.null) ) ])
# 	di_web_all=c(di_web_all, di_web[(!sapply(di_web,is.null) ) ])
# 	te_web_all=c(te_web_all, te_web[(!sapply(te_web,is.null) ) ])
# 	si_web_all=c(si_web_all, si_web[(!sapply(si_web,is.null) ) ])
# 	out1_all=c(out1_all, out1[1:nwebs])

# }

# #Take variables out of the lists to plot: 
# ncells=length(si_web_all)
# rDIT = data.frame(matrix(0, nrow=ncells, ncol =7) ) 
# ncnames = c("fwno","Biomass", "var_Biomass", "nspp", "shannon", "MI", "AI" )
# colnames(rDIT) = ncnames

# for (n in 1:ncells){
# 	rDIT$fwno[n] = n 
# 	rDIT$nspp[n] = out1_all[[n]]$spp_prms$nspp #Number of species

# 	tlast = dim(out1_all[[n]]$out)[1] - 1 #Length of time series
# 	rDIT$Biomass[n] = sum(out1_all[[n]]$out[tlast, 2:(rDIT$nspp[n]+1) ]) #Biomass at last time

# 	tbck = tlast*3/4 #Use a subset that excludes transient stage for variance
# 	rDIT$var_Biomass[n] = var( rowSums( out1_all[[n]]$out[ (tlast-tbck):tlast, 2:(rDIT$nspp[n]+1) ]) )

# 	#Shannon Diversity
# 	pi = out1_all[[n]]$out[tlast, 2:(rDIT$nspp[n]+1) ] / rDIT$Biomass[n]
# 	pi[pi <= 0 ] = NA
# 	rDIT$shannon[n] = - sum( pi*log(pi),na.rm=T )

# 	#Sum of Active Information
# 	aip = (pi*di_web_all[[n]]$ai_means)
# 	rDIT$AI[n] = sum(aip, na.rm=T)
# }

# lm_nspp = lm(rDIT$Biomass~rDIT$nspp)
# lm_H=lm(rDIT$Biomass~rDIT$shannon)
# lm_AI=lm(rDIT$Biomass~rDIT$AI)

# #Plots 
# ggplot (rDIT, aes(x = nspp, y = Biomass,color = "1" ) ) + geom_point () + 
# 	geom_point (aes(x = shannon, y =Biomass,color = "2"))+
# 	geom_point( aes (x = AI, y=Biomass,color = "3" ) )+
# 	scale_color_discrete(name ="", labels = c("# Species", "SDI","AI" ) )

