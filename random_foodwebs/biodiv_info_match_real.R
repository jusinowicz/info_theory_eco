#=============================================================================
# This code pairs with the output from biodiv_info_foodwebs. 
# It takes the result of the simulations and combines all Resource species
# into a single unit, which matches the "real" foodwebs more closely in terms
# of available data. 
#=============================================================================
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
#
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
#Within this loop, combine resource species and make the new data variables to be
#analyzed. 

ncells=length(aiE_web_all)
for (n in 1:ncells){
	k=1
	tlast = dim(out1_all[[n]]$out)[1] - 1 #Length of time series
	nRsp =  out1_all[[n]]$spp_prms$nRsp
	nCsp = out1_all[[n]]$spp_prms$nCsp
	nPsp = out1_all[[n]]$spp_prms$nPsp
	nspp = nRsp+nCsp+nPsp
	

	#Make new out1_all_new$out: 
	if( nRsp > 1){
		out1_all_new = NULL
		out1_all_tmp = out1_all[[n]]$out[,(nRsp+2):(nspp+1)]
		out1_all_tmp = cbind(out1_all[[n]]$out[,1], rowSums(out1_all[[n]]$out[,(2:(nRsp+1) )]),out1_all_tmp  )
		out1_all_new$out=out1_all_tmp

		#Make new out1_all_new$spp_prms
		out1_all_new$spp_prms = out1_all[[n]]$spp_prms
		out1_all_new$spp_prms$nRsp = 1 
		out1_all_new$spp_prms$nspp =  out1_all_new$spp_prms$nRsp+out1_all_new$spp_prms$nCsp+out1_all_new$spp_prms$nPsp
		out1_all_new$spp_prms$rR = sum(out1_all[[n]]$spp_prms$rR )
		out1_all_new$spp_prms$Ki = sum(out1_all[[n]]$spp_prms$Ki )
		out1_all_new$spp_prms$eFc = out1_all[[n]]$spp_prms$eFc[,1] 
		out1_all_new$spp_prms$cC = as.matrix(colSums(out1_all[[n]]$spp_prms$cC ))

		tlast = dim(out1_all_new$out)[1] - 1 #Length of time series
		nRsp =  out1_all_new$spp_prms$nRsp
		nCsp = out1_all_new$spp_prms$nCsp
		nPsp = out1_all_new$spp_prms$nPsp
		nspp = nRsp+nCsp+nPsp
	}else{out1_all_new = out1_all[[n]]}


	nt1 = tlast - 100
	nt2 = tlast

	# rweb1_all[n] = list(rutledge_web( spp_list=c(nRsp,nCsp,nPsp), pop_ts = out1_all[[n]]$out[nt1:nt2,2:(nspp+1)],
	# 	spp_prms = out1_all[[n]]$spp_prms, if_conditional = FALSE) )

	#=============================================================================
	aiE_web_all[n] = list( get_ais (  series1 = floor(out1_all_new$out[nt1:nt2,2:(nspp+1)]), 
		k=k, ensemble = TRUE)    )

	#=============================================================================
	MMI_web_all[n] = list( get_ais (  series1 = floor(out1_all_new$out[nt1:nt2,2:(nspp+1)]), 
		k=k, ensemble = TRUE)    )

	#Run these without any Resource species for comparison with the real data: 
		spp_prms1 = out1_all_new$spp_prms
		spp_prms1$rR = spp_prms1$rR[nRsp]
		spp_prms1$Ki = spp_prms1$Ki[nRsp]
		spp_prms1$cC = spp_prms1$cC[nRsp,]

		rweb1_all[n] = list(rutledge_web( spp_list=c(1,nCsp,nPsp), pop_ts = out1_all_new$out[nt1:nt2,(1+nRsp):(nspp+1)],
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
rDIT_sim_match = data.frame(matrix(0, nrow=ncells, ncol =12) ) 
ncnames = c("fwno","Biomass", "var_Biomass", "Snspp", "Fnspp", "shannon", "rS","rCE","rMI", "MI", 
	"AI","eq_state" )
colnames(rDIT_sim_match) = ncnames

for (n in 1:ncells){
	tlast1 = dim(out1_all[[n]]$out)[1] - 1 #Length of time series
	tlast2 = dim(aiE_web_all[[n]]$local)[1] - 1 #Length of time series


	rDIT_sim_match$fwno[n] = n 
	rDIT_sim_match$Snspp[n] = out1_all[[n]]$spp_prms$nspp #Starting number of species

	#####This needs to match the code above -- Are the Rsp being counted or not? 
	nRsp = out1_all[[n]]$spp_prms$nRsp
	rDIT_sim_match$Fnspp[n] = sum(out1_all[[n]]$out[tlast1,] > 0) - nRsp  #Final number of species

	rDIT_sim_match$Biomass[n] = sum(out1_all[[n]]$out[tlast1, 2:(rDIT_sim_match$Snspp[n]+1) ]) #Biomass at last time

	tbck = 1 #tlast*3/4 #Use a subset that excludes transient stage for variance
	rDIT_sim_match$var_Biomass[n] = var( rowSums( out1_all[[n]]$out[ (tlast1-tbck):tlast, 2:(rDIT_sim_match$Snspp[n]+1) ]) )

	#Shannon Diversity
	pi = out1_all[[n]]$out[tlast1, 2:(rDIT_sim_match$Snspp[n]+1) ] / rDIT_sim_match$Biomass[n]
	pi[pi <= 0 ] = NA
	rDIT_sim_match$shannon[n] = - sum( pi*log(pi),na.rm=T )

	#Rutledge Shannon Diversity, Conditional Entropy, and Mutual Information:  
	rDIT_sim_match$rS[n] = rweb1_all[[n]]$sD[tlast2]
	rDIT_sim_match$rCE[n] = rweb1_all[[n]]$ce2[tlast2]
	rDIT_sim_match$rMI[n] = rweb1_all[[n]]$mI_mean[tlast2]

	#Multiple Mutual Information
	rDIT_sim_match$MI[n] = MMI_web_all[[n]]$mean

	#Ensemble active information
	rDIT_sim_match$AI[n] = aiE_web_all[[n]]$mean

	#Determine whether this was a web in equilibrium (0) or not (1).
	eqtf = factor(levels=c(0,1))
	eq_test = test_eq( foodweb = out1_all[[n]], eqtest =tlast-50, t_type = "deriv")
	if(sum(eq_test)>1) {rDIT_sim_match$eq_state[n] = levels(eqtf)[2]} else { rDIT_sim_match$eq_state[n] = levels(eqtf)[1]   }


}

rDIT_sim_match_eq = subset(rDIT_sim_match, eq_state == 0)

#Log-log 
l_rDIT_sim_match_eq = log (rDIT_sim_match_eq[,1:9])
#l_rDIT_sim_match_eq[!is.finite(l_rDIT_sim_match_eq)] = NA
lm_nspp_sim_match = lm(l_rDIT_sim_match_eq$Biomass~l_rDIT_sim_match_eq$Fnspp)
lm_nspp_log_sim_match = lm(l_rDIT_sim_match_eq$Biomass~l_rDIT_sim_match_eq$Fnspp)
lm_H_sim_match=lm(l_rDIT_sim_match_eq$Biomass~l_rDIT_sim_match_eq$shannon)
lm_rS_sim_match=lm(l_rDIT_sim_match_eq$Biomass~l_rDIT_sim_match_eq$rS)
lm_rCE_sim_match=lm(l_rDIT_sim_match_eq$Biomass~l_rDIT_sim_match_eq$rCE)
lm_rMI_sim_match=lm(l_rDIT_sim_match_eq$Biomass~l_rDIT_sim_match_eq$rMI)
lm_rCEMI_sim_match=lm(l_rDIT_sim_match_eq$Biomass~l_rDIT_sim_match_eq$rCE+l_rDIT_sim_match_eq$rMI)
lm_rSMI_sim_match=lm(l_rDIT_sim_match_eq$Biomass~l_rDIT_sim_match_eq$rS+l_rDIT_sim_match_eq$rMI)

#Predict data from models to fit to figure

pr_nspp_log_sim_match = data.frame( Biomass =  exp(predict.lm ( lm_nspp_log_sim_match) ),
	Fnspp = exp(l_rDIT_sim_match_eq$Fnspp ) )
pr_H_sim_match = data.frame( Biomass = exp(predict( lm_H_sim_match) ) ,	shannon =exp(l_rDIT_sim_match_eq$shannon ) )
pr_rS_sim_match =data.frame( Biomass = exp(predict (lm_rS_sim_match ) ), rS= exp(l_rDIT_sim_match_eq$rS  ) )
pr_rCE_sim_match = data.frame( Biomass = exp(predict(lm_rCE_sim_match) ), rCE = exp(l_rDIT_sim_match_eq$rCE ) )
pr_rMI_sim_match = data.frame( Biomass = exp(predict(lm_rMI_sim_match) ), rMI = exp(l_rDIT_sim_match_eq$rMI ) )

#Plot of data with fitted lines: 
ggplot ( ) + 
	geom_point (data= rDIT_sim_match_eq, aes(x = (Fnspp), y = (Biomass),color = "1" )) + 
	geom_line ( data = pr_nspp_log_sim_match, aes(x = (Fnspp), y = (Biomass),color = "1" ) ) + 
	geom_point (data= rDIT_sim_match_eq,aes(x = (shannon), y = (Biomass),color = "2")) +
	geom_line ( data = pr_H_sim_match, aes(x = shannon, y = Biomass,color = "2" ) ) + 
	geom_point (data= rDIT_sim_match_eq,aes(x = (rS), y =(Biomass),color = "3")) +
	geom_line ( data = pr_rS_sim_match, aes(x = rS, y = Biomass,color = "3" ) ) + 
	geom_point( data= rDIT_sim_match_eq,aes (x = (rMI), y=(Biomass),color = "4" ) ) +
	geom_line ( data = pr_rMI_sim_match, aes(x = rMI, y = Biomass,color = "4" ) ) + 
	geom_point( data= rDIT_sim_match_eq,aes (x = (rCE), y=(Biomass),color = "5" ) ) +
	geom_line ( data = pr_rCE_sim_match, aes(x = rCE, y = Biomass,color = "5" ) ) + 
	scale_y_log10()+ scale_x_log10() +
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("# Species", "SDI", "rSD", "rMI","rCE") )
