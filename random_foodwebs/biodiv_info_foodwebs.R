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
source("../info_theory_functions/database_functions.R")



#=============================================================================
# Outer loop. Set the number of trials and determine how to generate 
# combinations of species and parameters. 
#=============================================================================

#Length and time steps of each model run
tend = 25
delta1 = 0.01
tl=tend/delta1

#The maximum block depth for dynamic info metrics (larger is more accurate, but
#slower and could cause crashing if too large)
k= 5 

#Number of food webs to generate
nwebs = 1

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
c1 = 1 
amp1 = 0.25 #100000 #1/exp(1) 

#Random consumers
c2 = 1
amp2 = 0.1
res_R = c(amp1,c1,amp2,c2)

for (w in 1:nwebs){ 
	print(w)

	#Assume 3 trophic levels unless otherwise specified.
	nRsp = ceiling(runif(1)*4)
	nCsp = ceiling(runif(1)*2)
	nPsp = ceiling(runif(1)*1)
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

		out1[w] = list(food_web_dynamics (spp_list = c(nRsp,nCsp,nPsp), spp_prms = spp_prms, 
		tend, delta1, res_R = NULL) )
		# print( paste( "nRsp", sum(out1[[w]]$out[tl,1:nRsp]>1) ) )
		# print( paste( "nCsp", sum(out1[[w]]$out[tl,(nRsp+1):nCsp]>1) ) )
		# print( paste( "nPsp", sum(out1[[w]]$out[tl,(nCsp+1):nPsp]>1) ) )		

		plot(out1[[w]]$out[,2], t="l", ylim = c(0, max(out1[[w]]$out[tl,],na.rm=T) ),col="red" )
		for(n in 2:(nRsp+1)) { lines(out1[[w]]$out[,n], col ="red") }
		for(n in (nRsp+2):(nRsp+nCsp+1) ){ lines(out1[[w]]$out[,n], col ="blue") }
		for(n in (nRsp+nCsp+1):(nspp+1) ){ lines(out1[[w]]$out[,n]) }

	
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
	#nt1 = 2/3*tl
	nt1 = tl - 100
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
#save(file = "rand_fwebmod7G.var", out1, rweb1,aiE_web,MMI_web) #These are deterministic
save(file = "/Volumes/TOSHIBA\ EXT/backups/mac_2020/Documents/GitHub/info_theory_eco/random_foodwebs/rand_fwebmod8C.var", out1, rweb1,aiE_web,MMI_web) #These are stochastic


#=============================================================================
# Load saved foodwebs and look at relationships between function and various 
# measures of complexity. 
#=============================================================================
#load("/Volumes/TOSHIBA\ EXT/backups/mac_2020/Documents/rand_fwebmod7G.var")
#This requires user input!
variable.list=list("out1", "di_web", "te_web","si_web")


file.name.list=c(
 # "/Volumes/TOSHIBA\ EXT/backups/mac_2020/Documents/GitHub/info_theory_eco/random_foodwebs/rand_fwebmod7A.var", 
 # "/Volumes/TOSHIBA\ EXT/backups/mac_2020/Documents/GitHub/info_theory_eco/random_foodwebs/rand_fwebmod7B.var",  
 # "/Volumes/TOSHIBA\ EXT/backups/mac_2020/Documents/GitHub/info_theory_eco/random_foodwebs/rand_fwebmod7C.var",
 # "/Volumes/TOSHIBA\ EXT/backups/mac_2020/Documents/GitHub/info_theory_eco/random_foodwebs/rand_fwebmod7D.var",
  "/Volumes/TOSHIBA\ EXT/backups/mac_2020/Documents/GitHub/info_theory_eco/random_foodwebs/rand_fwebmod7E.var",
  "/Volumes/TOSHIBA\ EXT/backups/mac_2020/Documents/GitHub/info_theory_eco/random_foodwebs/rand_fwebmod7F.var",
  "/Volumes/TOSHIBA\ EXT/backups/mac_2020/Documents/GitHub/info_theory_eco/random_foodwebs/rand_fwebmod7G.var"
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
	nwebs = (!sapply(aiE_web,is.null) )
	rweb1_all=c(rweb1_all, rweb1[nwebs ])
	aiE_web_all=c(aiE_web_all, aiE_web[nwebs])
	MMI_web_all=c(MMI_web_all, MMI_web[nwebs ])
	out1_all=c(out1_all, out1[nwebs])

}

#Re-run the loaded variables over a different range if need be: 
ncells=length(aiE_web_all)
for (n in 1:ncells){
	k=1
	tlast = dim(out1_all[[n]]$out)[1] - 1 #Length of time series
	nRsp = out1_all[[n]]$spp_prms$nRsp
	nCsp = out1_all[[n]]$spp_prms$nCsp
	nPsp = out1_all[[n]]$spp_prms$nPsp
	nspp = nRsp+nCsp+nPsp
	
	nt1 = tlast - 100
	nt2 = tlast

	# rweb1_all[n] = list(rutledge_web( spp_list=c(nRsp,nCsp,nPsp), pop_ts = out1_all[[n]]$out[nt1:nt2,2:(nspp+1)],
	# 	spp_prms = out1_all[[n]]$spp_prms, if_conditional = FALSE) )

	#=============================================================================
	aiE_web_all[n] = list( get_ais (  series1 = floor(out1_all[[n]]$out[nt1:nt2,2:(nspp+1)]), 
		k=k, ensemble = TRUE)    )

	#=============================================================================
	MMI_web_all[n] = list( get_ais (  series1 = floor(out1_all[[n]]$out[nt1:nt2,2:(nspp+1)]), 
		k=k, ensemble = TRUE)    )

	#Run these without any Resource species for comparison with the real data: 
		spp_prms1 = out1_all[[n]]$spp_prms
		spp_prms1$rR = spp_prms1$rR[nRsp]
		spp_prms1$Ki = spp_prms1$Ki[nRsp]
		spp_prms1$cC = spp_prms1$cC[nRsp,]

		rweb1_all[n] = list(rutledge_web( spp_list=c(1,nCsp,nPsp), pop_ts = out1_all[[n]]$out[nt1:nt2,(1+nRsp):(nspp+1)],
		spp_prms = spp_prms1, if_conditional = FALSE) )

	#=============================================================================
	# aiE_web_all[n] = list( get_ais (  series1 = floor(out1_all[[n]]$out[nt1:nt2,(1+nRsp):(nspp+1)]), 
	# 	k=k, ensemble = TRUE)    )

	# #=============================================================================
	# MMI_web_all[n] = list( get_ais (  series1 = floor(out1_all[[n]]$out[nt1:nt2,(1+nRsp):(nspp+1)]), 
	# 	k=k, ensemble = TRUE)    )
}

#Take variables out of the lists to plot: 
ncells=length(aiE_web_all)
rDIT_sim = data.frame(matrix(0, nrow=ncells, ncol =12) ) 
ncnames = c("fwno","Biomass", "var_Biomass", "Snspp", "Fnspp", "shannon", "rS","rCE","rMI", "MI", 
	"AI","eq_state" )
colnames(rDIT_sim) = ncnames

for (n in 1:ncells){
	tlast1 = dim(out1_all[[n]]$out)[1] - 1 #Length of time series
	tlast2 = dim(aiE_web_all[[n]]$local)[1] - 1 #Length of time series


	rDIT_sim$fwno[n] = n 
	rDIT_sim$Snspp[n] = out1_all[[n]]$spp_prms$nspp #Starting number of species

	#####This needs to match the code above -- Are the Rsp being counted or not? 
	nRsp = out1_all[[n]]$spp_prms$nRsp
	rDIT_sim$Fnspp[n] = sum(out1_all[[n]]$out[tlast1,] > 0) - nRsp-1  #Final number of species

	rDIT_sim$Biomass[n] = sum(out1_all[[n]]$out[tlast1, 2:(rDIT_sim$Snspp[n]+1) ]) #Biomass at last time

	tbck = 1 #tlast*3/4 #Use a subset that excludes transient stage for variance
	rDIT_sim$var_Biomass[n] = var( rowSums( out1_all[[n]]$out[ (tlast1-tbck):tlast, 2:(rDIT_sim$Snspp[n]+1) ]) )

	#Shannon Diversity
	pi = out1_all[[n]]$out[tlast1, 2:(rDIT_sim$Snspp[n]+1) ] / rDIT_sim$Biomass[n]
	pi[pi <= 0 ] = NA
	rDIT_sim$shannon[n] = - sum( pi*log(pi),na.rm=T )

	#Rutledge Shannon Diversity, Conditional Entropy, and Mutual Information:  
	rDIT_sim$rS[n] = rweb1_all[[n]]$sD[tlast2]
	rDIT_sim$rCE[n] = rweb1_all[[n]]$ce2[tlast2]
	rDIT_sim$rMI[n] = rweb1_all[[n]]$mI_mean[tlast2]

	#Multiple Mutual Information
	rDIT_sim$MI[n] = MMI_web_all[[n]]$mean

	#Ensemble active information
	rDIT_sim$AI[n] = aiE_web_all[[n]]$mean

	#Determine whether this was a web in equilibrium (0) or not (1).
	eqtf = factor(levels=c(0,1))
	eq_test = test_eq( foodweb = out1_all[[n]], eqtest =tlast-50, t_type = "deriv")
	if(sum(eq_test)>1) {rDIT_sim$eq_state[n] = levels(eqtf)[2]} else { rDIT_sim$eq_state[n] = levels(eqtf)[1]   }


}

rDIT_sim_eq = subset(rDIT_sim, eq_state == 0)

#Log-log 
l_rDIT_sim_eq = log (rDIT_sim_eq[,1:9])
#l_rDIT_sim_eq[!is.finite(l_rDIT_sim_eq)] = NA
lm_nspp_sim = lm(l_rDIT_sim_eq$Biomass~l_rDIT_sim_eq$Fnspp)
lm_nspp_log_sim = lm(l_rDIT_sim_eq$Biomass~l_rDIT_sim_eq$Fnspp)
lm_H_sim=lm(l_rDIT_sim_eq$Biomass~l_rDIT_sim_eq$shannon)
lm_rS_sim=lm(l_rDIT_sim_eq$Biomass~l_rDIT_sim_eq$rS)
lm_rCE_sim=lm(l_rDIT_sim_eq$Biomass~l_rDIT_sim_eq$rCE)
lm_rMI_sim=lm(l_rDIT_sim_eq$Biomass~l_rDIT_sim_eq$rMI)
lm_rCEMI_sim=lm(l_rDIT_sim_eq$Biomass~l_rDIT_sim_eq$rCE+l_rDIT_sim_eq$rMI)
lm_rSMI_sim=lm(l_rDIT_sim_eq$Biomass~l_rDIT_sim_eq$rS+l_rDIT_sim_eq$rMI)

#Predict data from models to fit to figure

pr_nspp_log_sim = data.frame( Biomass =  exp(predict.lm ( lm_nspp_log_sim) ),
	Fnspp = exp(l_rDIT_sim_eq$Fnspp ) )
pr_H_sim = data.frame( Biomass = exp(predict( lm_H_sim) ) ,	shannon =exp(l_rDIT_sim_eq$shannon ) )
pr_rS_sim =data.frame( Biomass = exp(predict (lm_rS_sim ) ), rS= exp(l_rDIT_sim_eq$rS  ) )
pr_rCE_sim = data.frame( Biomass = exp(predict(lm_rCE_sim) ), rCE = exp(l_rDIT_sim_eq$rCE ) )
pr_rMI_sim = data.frame( Biomass = exp(predict(lm_rMI_sim) ), rMI = exp(l_rDIT_sim_eq$rMI ) )

#Plot of data with fitted lines: 
ggplot ( ) + 
	geom_point (data= rDIT_sim_eq, aes(x = (Fnspp), y = (Biomass),color = "1" )) + 
	geom_line ( data = pr_nspp_log_sim, aes(x = (Fnspp), y = (Biomass),color = "1" ) ) + 
	geom_point (data= rDIT_sim_eq,aes(x = (shannon), y = (Biomass),color = "2")) +
	geom_line ( data = pr_H_sim, aes(x = shannon, y = Biomass,color = "2" ) ) + 
	geom_point (data= rDIT_sim_eq,aes(x = (rS), y =(Biomass),color = "3")) +
	geom_line ( data = pr_rS_sim, aes(x = rS, y = Biomass,color = "3" ) ) + 
	geom_point( data= rDIT_sim_eq,aes (x = (rMI), y=(Biomass),color = "4" ) ) +
	geom_line ( data = pr_rMI_sim, aes(x = rMI, y = Biomass,color = "4" ) ) + 
	geom_point( data= rDIT_sim_eq,aes (x = (rCE), y=(Biomass),color = "5" ) ) +
	geom_line ( data = pr_rCE_sim, aes(x = rCE, y = Biomass,color = "5" ) ) + 
	scale_y_log10()+ scale_x_log10() +
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("# Species", "SDI", "rSD", "rMI","rCE") )

ggsave("./rsVbio_rands1_sub.pdf", width = 8, height = 10)



#Log-y
lm_nspp_sim = lm(rDIT_sim_eq$Biomass~rDIT_sim_eq$Fnspp)
lm_nspp_log_sim = lm(I(log(rDIT_sim_eq$Biomass))~rDIT_sim_eq$Fnspp)
lm_H_sim=lm(I(log(rDIT_sim_eq$Biomass))~rDIT_sim_eq$shannon)
lm_rS_sim=lm(I(log(rDIT_sim_eq$Biomass))~rDIT_sim_eq$rS)
lm_rCE_sim=lm(I(log(rDIT_sim_eq$Biomass))~rDIT_sim_eq$rCE)
lm_rMI_sim=lm(I(log(rDIT_sim_eq$Biomass))~rDIT_sim_eq$rMI)
lm_rCEMI_sim=lm(I(log(rDIT_sim_eq$Biomass))~rDIT_sim_eq$rCE+rDIT_sim_eq$rMI)
lm_rSMI_sim=lm(I(log(rDIT_sim_eq$Biomass))~rDIT_sim_eq$rS+rDIT_sim_eq$rMI)

summary(lm_nspp_sim )
summary(lm_nspp_log_sim )
summary(lm_H_sim )
summary(lm_rS_sim )
summary(lm_rCE_sim)
summary(lm_rMI_sim )
summary(lm_rCEMI_sim )
summary(lm_rSMI_sim )

#Predict data from models to fit to figure
pr_nspp_sim = data.frame( Biomass = coef(lm_nspp_sim)[1] + coef(lm_nspp_sim)[2]* (1:max(rDIT_sim_eq$Fnspp) ),
 Fnspp = 1:max(rDIT_sim_eq$Fnspp) )
pr_nspp_log_sim = data.frame( Biomass = exp(coef(lm_nspp_log_sim)[1]) * exp(coef(lm_nspp_log_sim)[2]*(1:max(rDIT_sim_eq$Fnspp)) ),
 Fnspp = (1:max(rDIT_sim_eq$Fnspp) ) )
pr_H_sim = data.frame( Biomass = exp(coef(lm_H_sim)[1] ) * exp( coef(lm_H_sim)[2]* seq(0.1,5,0.1) ) ,
 shannon = seq(0.1,5,0.1) )
pr_rS_sim =data.frame( Biomass = exp(coef(lm_rS_sim)[1] ) * exp( coef(lm_rS_sim)[2]* seq(0.1,5,0.1) ) ,
 rS= seq(0.1,5,0.1) )
pr_rCE_sim = data.frame( Biomass = exp(coef(lm_rCE_sim)[1] ) * exp(coef(lm_rCE_sim)[2]* seq(0.1,5,0.1) ),
 rCE = seq(0.1,5,0.1) )
pr_rMI_sim = data.frame( Biomass = exp(coef(lm_rMI_sim)[1] ) * exp( coef(lm_rMI_sim)[2]* seq(0.1,5,0.1) ),
 rMI= seq(0.1,5,0.1) )
# pr_rCEMI
#pr_rSMI


#Predict data from models to fit to figure

# pr_nspp_log = data.frame( Biomass =  exp(predict.lm ( lm_nspp_log) ),
# 	Fnspp = exp(l_rDIT_sim_eq$Fnspp ) )
# pr_H = data.frame( Biomass = exp(predict( lm_H) ) ,	shannon =(l_rDIT_sim_eq$shannon ) )
# pr_rS =data.frame( Biomass = exp(predict (lm_rS ) ), rS= (l_rDIT_sim_eq$rS  ) )
# pr_rCE = data.frame( Biomass = exp(predict(lm_rCE) ), rCE = (l_rDIT_sim_eq$rCE ) )
# pr_rMI = data.frame( Biomass = exp(predict(lm_rMI) ), rMI = (l_rDIT_sim_eq$rMI ) )

#Plot of data with fitted lines: 
ggplot ( ) + 
	geom_point (data= rDIT_sim_eq, aes(x = (nspp), y = (Biomass),color = "1" )) + 
	geom_line ( data = pr_nspp_log_sim, aes(x = (nspp), y = (Biomass),color = "1" ) ) + 
	geom_point (data= rDIT_sim_eq,aes(x = (shannon), y = (Biomass),color = "2")) +
	geom_line ( data = pr_H_sim, aes(x = shannon, y = Biomass,color = "2" ) ) + 
	geom_point (data= rDIT_sim_eq,aes(x = (rS), y =(Biomass),color = "3")) +
	geom_line ( data = pr_rS_sim, aes(x = rS, y = Biomass,color = "3" ) ) + 
	geom_point( data= rDIT_sim_eq,aes (x = (rMI), y=(Biomass),color = "4" ) ) +
	geom_line ( data = pr_rMI_sim, aes(x = rMI, y = Biomass,color = "4" ) ) + 
	geom_point( data= rDIT_sim_eq,aes (x = (rCE), y=(Biomass),color = "5" ) ) +
	geom_line ( data = pr_rCE_sim, aes(x = rCE, y = Biomass,color = "5" ) ) + 
	scale_y_log10()+ 
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("# Species", "SDI", "rSD", "rMI","rCE") )



rDIT_sim_non = subset(rDIT_sim, eq_state == 0)

lm_Snspp = lm(rDIT_sim_eq$Biomass~rDIT_sim_eq$Snspp)
lm_Fnspp = lm(rDIT_sim_eq$Biomass~rDIT_sim_eq$Fnspp)
lm_H=lm(rDIT_sim_eq$Biomass~rDIT_sim_eq$shannon)
lm_rS=lm(I(log(rDIT_sim_eq$Biomass+1))~I(log(rDIT_sim_eq$rS+1) ))
lm_rCE=lm(I(log(rDIT_sim_eq$Biomass+1))~I(log(rDIT_sim_eq$rCE+1) ))
lm_rMI=lm(I(log(rDIT_sim_eq$Biomass+1))~I(log(rDIT_sim_eq$rMI+1) ))
lm_MI=lm(rDIT_sim_eq$Biomass~rDIT_sim_eq$MI)
lm_AI=lm(rDIT_sim_eq$Biomass~rDIT_sim_eq$AI)

summary(lm_Snspp )
summary(lm_Fnspp )
summary(lm_H )
summary(lm_rS )
summary(lm_rCE )
summary(lm_rMI )
summary(lm_MI )
summary(lm_AI )

lm_nspp2 = lm(rDIT_eq$Biomass~rDIT_eq$nspp+rDIT_eq$shannon)
lm_nspp3 = lm(rDIT_eq$Biomass~rDIT_eq$nspp+rDIT_eq$rMI)
lm_nspp4 = lm(rDIT_eq$Biomass~rDIT_eq$nspp+rDIT_eq$MI)


#Plots 
ggplot ( ) + 
	geom_point (data= rDIT_eq, aes(x = Fnspp, y = Biomass,color = "1" )) + 
	geom_point (data= rDIT_eq,aes(x = shannon, y =Biomass,color = "2")) +
	geom_point( data= rDIT_eq,aes (x = rMI, y=Biomass,color = "3" ) ) +
	geom_point( data= rDIT_eq,aes (x = MI, y=Biomass,color = "4" ) ) +
	geom_point( data= rDIT_eq,aes (x = AI, y=Biomass,color = "5" ) ) +
	scale_y_log10()+ scale_x_log10() +
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("# Species", "SDI", "rMI","MI","AI" ) )

#Plots 
ggplot ( ) + 
	#geom_point (data= rDIT_eq, aes(y = Fnspp, x = Biomass,color = "1" )) + 
	geom_point (data= rDIT_eq,aes(y = shannon, x =Biomass,color = "2")) +
	geom_point (data= rDIT_eq,aes(y = rS, x =Biomass,color = "3")) +
	geom_point (data= rDIT_eq,aes(y = rCE, x =Biomass,color = "4")) +
	geom_point( data= rDIT_eq,aes (y = rMI, x=Biomass,color = "5" ) ) +
	# geom_point( data= rDIT_eq,aes (x = MI, y=Biomass,color = "6" ) ) +
	# geom_point( data= rDIT_eq,aes (x = AI, y=Biomass,color = "7" ) ) +
	#scale_y_log10()+ scale_x_log10() +
	ylab("Complexity (Bits) ")+
	xlab("Biomass")+
	#scale_color_discrete(name ="", labels = c("# Species", "SDI", "rS","rCE", "rMI","MI","AI" ) )
	scale_color_discrete(name ="", labels = c("ShannonDI", "rShannon","rConditional Entropy",
	 "rMutual Information") )
ggsave("./complexity_v_biomass_all.pdf", width = 8, height = 10)

ggplot ( ) + 
	geom_point (data= rDIT, aes(x = nspp, y = var_Biomass,color = "1" )) + 
	geom_point (data= rDIT,aes(x = shannon, y =var_Biomass,color = "2")) +
	geom_point( data= rDIT,aes (x = rMI, y=var_Biomass,color = "3" ) ) +
	geom_point( data= rDIT,aes (x = MI, y=var_Biomass,color = "4" ) ) +
	geom_point( data= rDIT,aes (x = AI, y=var_Biomass,color = "5" ) ) +
	#scale_y_log10()+ scale_x_log10() +
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("# Species", "SDI", "rMI","MI","AI" ) )


lm_nsppv = lm(I(log(rDIT$var_Biomass))~I(log(rDIT$nspp)))
lm_Hv=lm(I(log(rDIT$var_Biomass))~I(log(rDIT$shannon)))
lm_rMIv=lm(I(log(rDIT$var_Biomass))~I(log(rDIT$rMI)))
lm_MIv=lm(I(log(rDIT$var_Biomass))~I(log(rDIT$MI)))
lm_AIv=lm(I(log(rDIT$var_Biomass))~I(log(rDIT$AI)))

summary(lm_nsppv )
summary(lm_Hv )
summary(lm_rMIv )
summary(lm_MIv )
summary(lm_AIv )

lm_nsppv2 = lm(I(log(rDIT$var_Biomass))~I(log(rDIT$nspp))+I(log(rDIT$shannon)))
lm_nsppv3 = lm(I(log(rDIT$var_Biomass))~I(log(rDIT$nspp))+I(log(rDIT$rMI)))
lm_nsppv4 = lm(I(log(rDIT$var_Biomass))~I(log(rDIT$nspp))+I(log(rDIT$MI)))

AIC(lm_nsppv,lm_nsppv2,lm_nsppv3,lm_nsppv4 )


#Check population plots: 
n=89
n2=61
n3=88
tlast = dim(out1_all[[n]]$out)[1] - 1
nRsp = out1_all[[n]]$spp_prms$nRsp
nCsp = out1_all[[n]]$spp_prms$nCsp
nPsp = out1_all[[n]]$spp_prms$nPsp

	plot(out1_all[[n]]$out[,1], t="l", ylim = c(0, max(out1_all[[n]]$out[tlast,],na.rm=T) ) )
		for(w in 2:nRsp){ lines(out1_all[[n]]$out[,w], col ="red") }
		for(w in (nRsp+2):(nRsp+nCsp+1) ){ lines(out1_all[[n]]$out[,w], col ="blue") }
		for(w in (nRsp+nCsp+2):(nRsp+nCsp+nPsp+1) ){ lines(out1_all[[n]]$out[,w]) }

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

