#=============================================================================
# R code to explore the Information Theoretic properties of real food webs. 
#	  
# 1. Download foodwebs from databases including EcoBase and Mangal. 
# 2. Use information theory to analyze the food-web complexity. 
# 3. Compare the explanatory power of IT complexity to Shannon Diversity and 
#    species number for biomass (ecosystem function)
#=============================================================================
#=============================================================================
# load libraries
#=============================================================================
library(tidyverse)
library(lubridate)
library(mgcv)
library(RCurl)
library(XML)
source("../info_theory_functions/food_web_functions.R")
source("../info_theory_functions/info_theory_functions.R")
source("../info_theory_functions/database_functions.R")

#=============================================================================
# EcoBase: Get all of the available full foodwebs
#=============================================================================
fwlist = get_eb() 
nwebs = length (fwlist)

#The maximum block depth for dynamic info metrics (larger is more accurate, but
#slower and could cause crashing if too large)
k= 5 

#Converting the web to Rutledge's compartment model and calculating the information
#theoretic quantities: Shannon Entropy, Mutual Information, Conditional Entropy
rweb_mb = vector("list",nwebs)

for (w in 1:nwebs){ 
	print(w)

	#Reformat these to a vector with named entries
	biomass = abs(fwlist[[w]]$biomass$biomass)
	names(biomass) = fwlist[[w]]$biomass$name 

	#Ratio of production P/B
	pb = fwlist[[w]]$pb$pb
	names(pb) = fwlist[[w]]$pb$name 

	#Ratio of consumption Q/B
	qb = fwlist[[w]]$qb$qb
	names(qb) = fwlist[[w]]$qb$name 

	#Ecotrophic efficiency EE
	ee = fwlist[[w]]$ee$ee
	names(ee) = fwlist[[w]]$ee$name 

	#The dietary matrix
	DC = fwlist[[w]]$trophic_relations

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
	
	rweb_mb[w] = list(rutledge_web_mb(biomass = biomass, pb = pb, qb = qb, DC = DC, ee = ee,  
						 if_conditional = FALSE) )

	
}

#save(file = "rand_fwebmod6F.var", out1,  di_web,te_web,si_web)
#save(file = "rand_fwebmod7G.var", out1, rweb1,aiE_web,MMI_web) #These are deterministic
#save(file = "/Volumes/TOSHIBA\ EXT/backups/mac_2020/Documents/GitHub/info_theory_eco/random_foodwebs/rand_fwebmod8C.var", rweb_mb) #These are stochastic


#Take variables out of the lists to plot: 
ncells=length(rweb_mb)
rDIT = data.frame(matrix(0, nrow=ncells, ncol =9) ) 
ncnames = c("fwno","Biomass", "var_Biomass", "nspp", "shannon", "rS","rCE","rMI")
colnames(rDIT) = ncnames

for (n in 1:ncells){
	rDIT$fwno[n] = n 
	rDIT$nspp[n] = length(rweb_mb[[n]]$biomass) #Final number of species

	rDIT$Biomass[n] = sum(rweb_mb[[n]]$biomass) #Total biomass

	rDIT$var_Biomass[n] = 0

	#Shannon Diversity
	pi = rweb_mb[[n]]$biomass / rDIT$Biomass[n]
	pi[pi <= 0 ] = NA
	rDIT$shannon[n] = - sum( pi*log(pi),na.rm=T )

	#Rutledge Shannon Diversity, Conditional Entropy, and Mutual Information:  
	rDIT$rS[n] = rweb_mb[[n]]$sD
	rDIT$rCE[n] = rweb_mb[[n]]$ce2
	rDIT$rMI[n] = rweb_mb[[n]]$mI_mean2

}

rDIT_eq = subset(rDIT, rMI >0.01 )

lm_nspp = lm(rDIT_eq$Biomass~rDIT_eq$nspp)
lm_H=lm(rDIT_eq$Biomass~rDIT_eq$shannon)
lm_rS=lm(I(log(rDIT_eq$Biomass))~I(log(rDIT_eq$rS) ))
lm_rCE=lm(I(log(rDIT_eq$Biomass))~I(log(rDIT_eq$rCE) ))
lm_rMI=lm(I(log(rDIT_eq$Biomass))~I(log(rDIT_eq$rMI) ))
lm_rCEMI=lm(I(log(rDIT_eq$Biomass))~I(log(rDIT_eq$rCE) )+I(log(rDIT_eq$rMI) ) )
lm_rSMI=lm(I(log(rDIT_eq$Biomass))~I(log(rDIT_eq$rS) )+I(log(rDIT_eq$rMI) ) )



summary(lm_nspp )
summary(lm_H )
summary(lm_rS )
summary(lm_rCE )
summary(lm_rMI )
summary(lm_rCEMI )
summary(lm_rSMI )


lm_nspp2 = lm(rDIT_eq$Biomass~rDIT_eq$nspp+rDIT_eq$shannon)
lm_nspp3 = lm(rDIT_eq$Biomass~rDIT_eq$nspp+rDIT_eq$rMI)
lm_nspp4 = lm(rDIT_eq$Biomass~rDIT_eq$nspp+rDIT_eq$MI)


#Plots 
ggplot ( ) + 
	geom_point (data= rDIT_eq, aes(x = nspp, y = Biomass,color = "1" )) + 
	geom_point (data= rDIT_eq,aes(x = shannon, y =Biomass,color = "2")) +
	geom_point (data= rDIT_eq,aes(x = rS, y =Biomass,color = "3")) +
	geom_point( data= rDIT_eq,aes (x = rMI, y=Biomass,color = "4" ) ) +
	geom_point( data= rDIT_eq,aes (x = rCE, y=Biomass,color = "5" ) ) +
	scale_y_log10()+ scale_x_log10() +
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("# Species", "SDI", "rSD", "rMI","rCE") )

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

