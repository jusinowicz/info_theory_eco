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
# Get all of the available full foodwebs. These are in CSVs with 2 prefixes: 
# Biomass has the biomass, P/B, Q/B, and EE
# DC has the dietary composition matrix
#=============================================================================
flist_b = list.files("./fwfw", pattern="iomass")  #Load the biomass 
flist_DC = list.files("./fwfw", pattern="DC")  #Load the biomass 
fwlist = get_eb_fw (biomass_file_list=flist_b, DC_file_list=flist_DC)

nwebs = length (fwlist)

#The maximum block depth for dynamic info metrics (larger is more accurate, but
#slower and could cause crashing if too large)
k= 5 

#Converting the web to Rutledge's compartment model and calculating the information
#theoretic quantities: Shannon Entropy, Mutual Information, Conditional Entropy
rweb_mb = vector("list",nwebs)
all_names = NULL

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
	ee[ee>1] = 1 #This should not happen, but it does. wtf???

	#The dietary matrix
	DC = fwlist[[w]]$trophic_relations

	#Percentage of consumption that is lost to Detritus
	gs = fwlist[[w]]$gs$gs
	names(gs) = fwlist[[w]]$gs$name

	#The balance is given by  pb*biomass*ee- rowSums(matrix(qb*biomass, nspp,nspp,byrow=T)*DC)

	#=======================
	#Do some pre-formatting for Detritus!!! Use the GS and Export of Detritus to 
	#calculate its inputs and outputs:
	#=======================

	#Total Detritus consumed: 
	nspp = length ( biomass )

	#Test whether there are multiple sub-detritus groups, as well as categorized as 
	# "discard" and "detached" 
	d_groups = biomass[grepl("etrit", names(biomass)) | grepl("iscard", names(biomass)) 
						 | grepl("etached", names(biomass)) | grepl("Dead", names(biomass)) 
						 | grepl("Matter", names(biomass)) ]
 	
 	n_d = length(d_groups)
	#Total detritus consumed by other organisms
	d_con = rowSums(matrix(qb*biomass, nspp,nspp,byrow=T)*DC)[
			names(biomass)%in%names(d_groups)]
		
 	#If they do not have a biomass value assigned to them (usually specified with -9999),
 	#then group all of the categories into one "Detritus" category and use the expected
 	#production of biomass from the GS of other organisms to assign it: 
 	#Test for the -9999 condition, which means biomass hasn't been specified
	test1 = abs(biomass[(names(biomass)%in%names(d_groups) )]) == 9999
	test2 = biomass[(names(biomass)%in%names(d_groups) )]< 0
	
	if( sum(test1) > 0 | sum(test2) > 0) {  
		
		biomass [names(biomass)%in%names(d_groups)] = d_con
	}
		#biomass [names(biomass)%in%names(d_groups)] = d_con
		

		# #Use gs as the Detritus "consumption" rates
		# DC[,"Detritus"] = gs/ sum (gs)
		# #QB:This produces a vector but every value should be the same
		# qb["Detritus"] = mean(gs/(DC[,"Detritus"] * biomass["Detritus"] ),na.rm=T)

		#Alternatively, make Detritus more like primary productivity: 
		DC[,names(d_groups)] = matrix(0,nspp,n_d)
		qb[names(d_groups)] = 0

		#PB 
		pb[names(d_groups)] = d_con/biomass[names(d_groups)] 
	
	
	
			#Calculate the replacement biomass field:
		# #Total Detritus produced by inefficency in consumption (GS) and mortality (1-EE)
		# d_prod = sum(biomass[(names(biomass)!="Detritus")]*
		# 	qb[(names(qb)!="Detritus")]*
		# 	gs[(names(gs)!="Detritus")]) + 
		# 	sum( (1-ee[(names(ee)!="Detritus")])*
		# 	biomass[(names(biomass)!="Detritus")])

		# d_prod = (biomass[!(names(biomass)%in%names(d_groups) )]*
		# 	qb[!(names(qb)%in%names(d_groups) )]*
		# 	gs[!(names(gs)%in%names(d_groups) )]) + 
		# 	( (1-ee[!(names(ee)%in%names(d_groups) )])*
		# 		biomass[!(names(biomass)%in%names(d_groups) )])

	#Always use d_prod for standing detritus biomass? 
	#biomass["Detritus"] = d_prod	

	
	#Remove Detritus: 
	# biomass= biomass[!(names(biomass)%in% names(d_groups) ) ]
	# qb= qb[!(names(qb)%in% names(d_groups) ) ]
	# pb= pb[!(names(pb)%in% names(d_groups) ) ]
	# ee= ee[!(names(ee)%in% names(d_groups) ) ]
	# DC = DC[!(colnames(DC)%in%names(d_groups) ),!(rownames(DC)%in%names(d_groups) )]


	#Remove Benthic categories:
	# b_groups = biomass[grepl("enthos", names(biomass)) | grepl("enthic", names(biomass)) ]
 	
 # 	n_b = length(b_groups)  
	# biomass= biomass[!(names(biomass)%in% names(b_groups) ) ]
	# qb= qb[!(names(qb)%in% names(b_groups) ) ]
	# pb= pb[!(names(pb)%in% names(b_groups) ) ]
	# ee= ee[!(names(ee)%in% names(b_groups) ) ]
	# DC = DC[!(colnames(DC)%in%names(b_groups) ),!(rownames(DC)%in%names(b_groups) )]


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

	all_names = c(all_names, names(biomass))

	
}

#save(file = "rand_fwebmod6F.var", out1,  di_web,te_web,si_web)
#save(file = "rand_fwebmod7G.var", out1, rweb1,aiE_web,MMI_web) #These are deterministic
#save(file = "/Volumes/TOSHIBA\ EXT/backups/mac_2020/Documents/GitHub/info_theory_eco/random_foodwebs/rand_fwebmod8C.var", rweb_mb) #These are stochastic


#Take variables out of the lists to plot: 
ncells=length(rweb_mb)
rDIT_lake = data.frame(matrix(0, nrow=ncells, ncol =9) ) 
ncnames = c("fwno","Biomass", "var_Biomass", "nspp", "shannon", "rS","rCE","rMI")
colnames(rDIT_lake) = ncnames

for (n in 1:ncells){
	rDIT_lake$fwno[n] = n 

	#Exclude Detritus from biomass count: 
	# bm1 = rweb_mb[[n]]$biomass[
	# 	names(rweb_mb[[n]]$biomass) != "Detritus" & 
	# 	names(rweb_mb[[n]]$biomass) != "detritus" ]
	
	#Don't exclude Detritus from biomass count: 
	bm1 = rweb_mb[[n]]$biomass
	
	rDIT_lake$nspp[n] = length(bm1) #Final number of species

	rDIT_lake$Biomass[n] = sum(bm1) #Total biomass

	rDIT_lake$var_Biomass[n] = 0

	#Shannon Diversity
	pi = bm1/sum(bm1)
	pi[pi <= 0 ] = NA
	rDIT_lake$shannon[n] = - sum( pi*log(pi),na.rm=T )

	#Rutledge Shannon Diversity, Conditional Entropy, and Mutual Information:  
	rDIT_lake$rS[n] = rweb_mb[[n]]$sD
	rDIT_lake$rCE[n] = rweb_mb[[n]]$ce2
	rDIT_lake$rMI[n] = rweb_mb[[n]]$mI_mean2

}


#Fit models to see what predicts Biomass the best
#rDIT_lake_eq = subset(rDIT_lake, rMI >0.01 )
rDIT_lake_eq = rDIT_lake #rDIT_lake[c(-5),]
rDIT_lake_eq = subset(rDIT_lake_eq, rMI >0.01 )

# rDIT_lake_eq_con = rDIT_lake_con[c(-71,-68),]
# rDIT_lake_eq_con = subset(rDIT_lake_eq_con, rMI >0.01 )

#Log-log 
l_rDIT_lake_eq = log (rDIT_lake_eq)
lm_nspp = lm(l_rDIT_lake_eq$Biomass~l_rDIT_lake_eq$nspp)
lm_nspp_log = lm(l_rDIT_lake_eq$Biomass~l_rDIT_lake_eq$nspp)
lm_H=lm(l_rDIT_lake_eq$Biomass~l_rDIT_lake_eq$shannon)
lm_rS=lm(l_rDIT_lake_eq$Biomass~l_rDIT_lake_eq$rS)
#lm_rS=lm(l_rDIT_lake_eq$Biomass~l_rDIT_lake_eq$rS+I(l_rDIT_lake_eq$rS^2))
lm_rCE=lm(l_rDIT_lake_eq$Biomass~l_rDIT_lake_eq$rCE)
lm_rMI=lm(l_rDIT_lake_eq$Biomass~l_rDIT_lake_eq$rMI)
lm_rCEMI=lm(l_rDIT_lake_eq$Biomass~l_rDIT_lake_eq$rCE+l_rDIT_lake_eq$rMI)
lm_rSMI=lm(l_rDIT_lake_eq$Biomass~l_rDIT_lake_eq$rS+l_rDIT_lake_eq$rMI)

gam_rS= gam (Biomass~s(rS),data=rDIT_lake_eq )
pr_gam_rS = predict(gam_rS, newdata =data.frame(rS = rDIT_lake_eq$rS), type="response")
#Predict data from models to fit to figure

pr_nspp_log = data.frame( Biomass =  exp(predict.lm ( lm_nspp_log) ),
	nspp = exp(l_rDIT_lake_eq$nspp ) )
pr_H = data.frame( Biomass = exp(predict( lm_H) ) ,	shannon =exp(l_rDIT_lake_eq$shannon ) )
pr_rS =data.frame( Biomass = exp(predict (lm_rS ) ), rS= exp(l_rDIT_lake_eq$rS  ) )
pr_rCE = data.frame( Biomass = exp(predict(lm_rCE) ), rCE = exp(l_rDIT_lake_eq$rCE ) )
pr_rMI = data.frame( Biomass = exp(predict(lm_rMI) ), rMI = exp(l_rDIT_lake_eq$rMI ) )


summary(lm_nspp )
summary(lm_nspp_log )
summary(lm_H )
summary(lm_rS )
summary(lm_rCE )
summary(lm_rMI )
summary(lm_rCEMI )
summary(lm_rSMI )

#Plots 

#Plot of data on log scale: 
ggplot ( ) + 
	geom_point (data= rDIT_lake_eq, aes(x = nspp, y = Biomass,color = "1" )) + 
	geom_point (data= rDIT_lake_eq,aes(x = shannon, y =Biomass,color = "2")) +
	geom_point (data= rDIT_lake_eq,aes(x = rS, y =Biomass,color = "3")) +
	geom_point( data= rDIT_lake_eq,aes (x = rMI, y=Biomass,color = "4" ) ) +
	geom_point( data= rDIT_lake_eq,aes (x = rCE, y=Biomass,color = "5" ) ) +
	scale_y_log10()+ scale_x_log10() +
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("# Species", "SDI", "rSD", "rMI","rCE") )

#Plot of data with fitted lines: 
ggplot ( ) + 
	geom_point (data= rDIT_lake_eq, aes(x = (nspp), y = (Biomass),color = "1" )) + 
	geom_line ( data = pr_nspp_log, aes(x = (nspp), y = (Biomass),color = "1" ) ) + 
	geom_point (data= rDIT_lake_eq,aes(x = (shannon), y = (Biomass),color = "2")) +
	geom_line ( data = pr_H, aes(x = shannon, y = Biomass,color = "2" ) ) + 
	geom_point (data= rDIT_lake_eq,aes(x = (rS), y =(Biomass),color = "3")) +
	geom_line ( data = pr_rS, aes(x = rS, y = Biomass,color = "3" ) ) + 
	geom_point( data= rDIT_lake_eq,aes (x = (rMI), y=(Biomass),color = "4" ) ) +
	geom_line ( data = pr_rMI, aes(x = rMI, y = Biomass,color = "4" ) ) + 
	geom_point( data= rDIT_lake_eq,aes (x = (rCE), y=(Biomass),color = "5" ) ) +
	geom_line ( data = pr_rCE, aes(x = rCE, y = Biomass,color = "5" ) ) + 
	scale_y_log10()+ scale_x_log10() +
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("# Species", "SDI", "rSD", "rMI","rCE") )

##================================
# Joint plots of real and sim data with fitted lines: 
##================================
#Species
ggplot ( ) + 
	geom_point (data= rDIT_lake_eq, aes(x = (nspp), y = (Biomass),color = "1" )) + 
	geom_text(data= rDIT_lake_eq, aes(x = (nspp), y = (Biomass),label= fwno),hjust=0, vjust=0)+
	geom_line ( data = pr_nspp_log, aes(x = (nspp), y = (Biomass),color = "1" ) ) + 
	geom_point (data= rDIT_lake_sim_eq, aes(x = (Fnspp), y = (Biomass),color = "2" )) + 
	geom_line ( data = pr_nspp_log_sim, aes(x = (Fnspp), y = (Biomass),color = "2" ) ) + 
	scale_y_log10()+ scale_x_log10() +
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("# Species", "#Species Sims") )

ggsave("./nsppVbio_rands1_sub.pdf", width = 8, height = 10)

#Shannon Diversity
ggplot ( ) + 
	geom_point (data= rDIT_lake_eq,aes(x = (shannon), y = (Biomass),color = "2")) +
	geom_line ( data = pr_H, aes(x = shannon, y = Biomass,color = "2" ) ) + 
	geom_point (data= rDIT_lake_sim_eq,aes(x = (shannon), y = (Biomass),color = "3")) +
	geom_line ( data = pr_H_sim, aes(x = shannon, y = Biomass,color = "3" ) ) + 
	scale_y_log10()+ scale_x_log10() +
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("SDI", "SDI Sims") )

ggsave("./sdiVbio_rands1.pdf", width = 8, height = 10)

#Rutledge Shannon
ggplot ( ) + 
	geom_point (data= rDIT_lake_eq,aes(x = (rS), y =(Biomass),color = "4")) +
	geom_text(data= rDIT_lake_eq, aes(x = (rS), y = (Biomass),label= fwno),hjust=0, vjust=0)+
	geom_line ( data = pr_rS, aes(x = rS, y = Biomass,color = "4" ) ) + 
	#geom_line ( data = pr_gam_rS, aes(x = rS, y = Biomass,color = "4" ) ) + 
	geom_point (data= rDIT_lake_sim_eq,aes(x = (rS), y =(Biomass),color = "5")) +
	geom_line ( data = pr_rS_sim, aes(x = rS, y = Biomass,color = "5" ) ) + 
	scale_y_log10()+ scale_x_log10() + 
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("rMI", "rMI Sims") )

ggsave("./rsVbio_rands1_sub.pdf", width = 8, height = 10)

#MI
ggplot ( ) + 
	geom_point( data= rDIT_lake_eq,aes (x = (rMI), y=(Biomass),color = "1" ) ) +
	geom_text(data= rDIT_lake_eq, aes(x = (rMI), y = (Biomass),label= fwno),hjust=0, vjust=0)+
	geom_line ( data = pr_rMI, aes(x = rMI, y = Biomass,color = "1" ) ) + 
	geom_point( data= rDIT_lake_sim_eq,aes (x = (rMI), y=(Biomass),color = "2" ) ) +
	geom_line ( data = pr_rMI_sim, aes(x = rMI, y = Biomass,color = "2" ) ) + 
	scale_y_log10()+ scale_x_log10() +
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("rMI", "rMI Sims") )

ggsave("./rMIVbio_rands1_sub.pdf", width = 8, height = 10)


#CE
ggplot()+ 
	geom_point( data= rDIT_lake_eq,aes (x = (rCE), y=(Biomass),color = "8" ) ) +
	geom_text(data= rDIT_lake_eq, aes(x = (rCE), y = (Biomass),label= fwno),hjust=0, vjust=0)+
	geom_line ( data = pr_rCE, aes(x = rCE, y = Biomass,color = "8" ) ) + 
	geom_point( data= rDIT_lake_sim_eq,aes (x = (rCE), y=(Biomass),color = "9" ) ) +
	geom_line ( data = pr_rCE_sim, aes(x = rCE, y = Biomass,color = "9" ) ) + 
	scale_y_log10()+ scale_x_log10() +
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("rCE", "rCE Sims") )

ggsave("./rCEVbio_rands1_nod.pdf", width = 8, height = 10)

####Distribution of biomass across species: 
# rDIT_lake_eq = rDIT_lake[c(-71,-68),]
# rDIT_lake_eq = subset(rDIT_lake_eq, rMI >0.01 )
l_rDIT_lake_eq = log (rDIT_lake_eq)

lm_sim_Bnspp = lm( I(log(rDIT_lake_sim_eq$Biomass/rDIT_lake_sim_eq$Fnspp) )~log(rDIT_lake_sim_eq$Fnspp) )
pr_sim_Bnspp = data.frame( bns_ratio = exp(predict( lm_sim_Bnspp) ) ,	Fnspp=exp(l_rDIT_lake_sim_eq$Fnspp) )

lm_Bnspp = lm( I(log(rDIT_lake_eq$Biomass/rDIT_lake_eq$nspp))~log(rDIT_lake_eq$nspp ) )
pr_Bnspp = data.frame( bns_ratio = exp(predict( lm_Bnspp) ) ,	nspp=exp(l_rDIT_lake_eq$nspp) )


ggplot()+ 
	geom_histogram( data= rDIT_lake_eq,aes (x=(Biomass)/nspp, color = "8" ) ) +
	geom_histogram( data= rDIT_lake_sim_eq,aes ( x=(Biomass)/Fnspp,color = "9" ) ) 

ggplot()+ 
	geom_point( data= rDIT_lake_eq,aes (x = nspp, y=(Biomass)/nspp, color = "8" ) ) +
	geom_line ( data = pr_Bnspp, aes(x = nspp, y = bns_ratio,color = "8" ) ) + 
	#geom_text( data= rDIT_lake_eq, aes(x = nspp, y = (Biomass)/nspp,label= fwno),hjust=0, vjust=0)+
	geom_point( data= (rDIT_lake_sim_eq),aes (x = Fnspp, y=(Biomass)/Fnspp,color = "9" ) ) +
	geom_line ( data = pr_sim_Bnspp, aes(x = Fnspp, y = bns_ratio,color = "9" ) ) + 
	scale_y_log10()+ scale_x_log10() +
	#ylim(0,100)+
	xlab("N Species")+
	ylab("Biomass/species")+
	scale_color_discrete(name ="", labels = c("real", "sims") )

ggsave("./biomass_per_species_sub.pdf", width = 8, height = 10)


lm_sim_Brce = lm( I(log(rDIT_lake_sim_eq$Biomass/rDIT_lake_sim_eq$rCE) )~log(rDIT_lake_sim_eq$rCE) )
pr_sim_Brce = data.frame( bns_ratio = exp(predict( lm_sim_Brce) ) ,	Fnspp=exp(l_rDIT_lake_sim_eq$rCE) )

lm_Brce = lm( I(log(rDIT_lake_eq$Biomass/rDIT_lake_eq$rCE))~log(rDIT_lake_eq$rCE ) )
pr_Brce = data.frame( bns_ratio = exp(predict( lm_Brce) ) ,	nspp=exp(l_rDIT_lake_eq$rCE) )


ggplot()+ 
	geom_point( data= rDIT_lake_eq,aes (x = rCE, y=(Biomass)/rCE, color = "8" ) ) +
	geom_line ( data = pr_Brce, aes(x = nspp, y = bns_ratio,color = "8" ) ) + 
	#geom_text( data= rDIT_lake_eq, aes(x = nspp, y = (Biomass)/nspp,label= fwno),hjust=0, vjust=0)+
	geom_point( data= (rDIT_lake_sim_eq),aes (x = rCE, y=(Biomass)/rCE,color = "9" ) ) +
	geom_line ( data = pr_sim_Brce, aes(x = Fnspp, y = bns_ratio,color = "9" ) ) + 
	scale_y_log10()+ scale_x_log10() +
	#ylim(0,100)+
	xlab("N Species")+
	ylab("Biomass/species")+
	scale_color_discrete(name ="", labels = c("real", "sims") )

ggsave("./biomass_per_CE_sub.pdf", width = 8, height = 10)

#Multiple
#Plot of data with fitted lines: 
ggplot ( ) + 
	geom_point (data= rDIT_lake_eq,aes(x = (shannon), y = (Biomass),color = "2")) +
	geom_line ( data = pr_H, aes(x = shannon, y = Biomass,color = "2" ) ) + 
	geom_point (data= rDIT_lake_sim_eq,aes(x = (shannon), y = (Biomass),color = "3")) +
	geom_line ( data = pr_H_sim, aes(x = shannon, y = Biomass,color = "3" ) ) + 
	geom_point (data= rDIT_lake_eq,aes(x = (rS), y =(Biomass),color = "4")) +
	geom_line ( data = pr_rS, aes(x = rS, y = Biomass,color = "4" ) ) + 
	geom_point (data= rDIT_lake_sim_eq,aes(x = (rS), y =(Biomass),color = "5")) +
	geom_line ( data = pr_rS_sim, aes(x = rS, y = Biomass,color = "5" ) ) + 
	geom_point( data= rDIT_lake_eq,aes (x = (rMI), y=(Biomass),color = "6" ) ) +
	geom_line ( data = pr_rMI, aes(x = rMI, y = Biomass,color = "6" ) ) + 
	geom_point( data= rDIT_lake_sim_eq,aes (x = (rMI), y=(Biomass),color = "7" ) ) +
	geom_line ( data = pr_rMI_sim, aes(x = rMI, y = Biomass,color = "7" ) ) + 
	geom_point( data= rDIT_lake_eq,aes (x = (rCE), y=(Biomass),color = "8" ) ) +
	geom_line ( data = pr_rCE, aes(x = rCE, y = Biomass,color = "8" ) ) + 
	geom_point( data= rDIT_lake_sim_eq,aes (x = (rCE), y=(Biomass),color = "9" ) ) +
	geom_line ( data = pr_rCE_sim, aes(x = rCE, y = Biomass,color = "9" ) ) + 
	scale_y_log10()+ scale_x_log10() +
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("SDI", "SDI Sim",
		"rSD", "rSD Sim","rMI","rMI Sim","rCE","rCE Sim") )


#Log-y
lm_nspp = lm(rDIT_lake_eq$Biomass~rDIT_lake_eq$nspp)
lm_nspp_log = lm(I(log(rDIT_lake_eq$Biomass))~rDIT_lake_eq$nspp)
lm_H=lm(I(log(rDIT_lake_eq$Biomass))~rDIT_lake_eq$shannon)
lm_rS=lm(I(log(rDIT_lake_eq$Biomass))~rDIT_lake_eq$rS)
lm_rCE=lm(I(log(rDIT_lake_eq$Biomass))~rDIT_lake_eq$rCE)
lm_rMI=lm(I(log(rDIT_lake_eq$Biomass))~rDIT_lake_eq$rMI)
lm_rCEMI=lm(I(log(rDIT_lake_eq$Biomass))~rDIT_lake_eq$rCE+rDIT_lake_eq$rMI)
lm_rSMI=lm(I(log(rDIT_lake_eq$Biomass))~rDIT_lake_eq$rS+rDIT_lake_eq$rMI)

summary(lm_nspp )
summary(lm_nspp_log )
summary(lm_H )
summary(lm_rS )
summary(lm_rCE )
summary(lm_rMI )
summary(lm_rCEMI )
summary(lm_rSMI )

gam_rS= gam (I(log(Biomass))~s(rS),data=rDIT_lake_eq )
pr_gam_rS = predict(gam_rS, newdata =data.frame(rS = rDIT_lake_eq$rS), type="response")
#Predict data from models to fit to figure
pr_nspp = data.frame( Biomass = coef(lm_nspp)[1] + coef(lm_nspp)[2]* (1:max(rDIT_lake_eq$nspp) ),
 nspp = 1:max(rDIT_lake_eq$nspp) )
pr_nspp_log = data.frame( Biomass = exp(coef(lm_nspp_log)[1]) * exp(coef(lm_nspp_log)[2]*(1:max(rDIT_lake_eq$nspp)) ),
 nspp = (1:max(rDIT_lake_eq$nspp) ) )
pr_H = data.frame( Biomass = exp(coef(lm_H)[1] ) * exp( coef(lm_H)[2]* seq(0.1,5,0.1) ) ,
 shannon = seq(0.1,5,0.1) )
pr_rS =data.frame( Biomass = exp(coef(lm_rS)[1] ) * exp( coef(lm_rS)[2]* seq(0.1,5,0.1) ) ,
 rS= seq(0.1,5,0.1) )
pr_rCE = data.frame( Biomass = exp(coef(lm_rCE)[1] ) * exp(coef(lm_rCE)[2]* seq(0.1,5,0.1) ),
 rCE = seq(0.1,5,0.1) )
pr_rMI = data.frame( Biomass = exp(coef(lm_rMI)[1] ) * exp( coef(lm_rMI)[2]* seq(0.1,5,0.1) ),
 rMI= seq(0.1,5,0.1) )
# pr_rCEMI
# pr_rSMI

#Plot of data with fitted lines: 
ggplot ( ) + 
	geom_point (data= rDIT_lake_eq, aes(x = (nspp), y = (Biomass),color = "1" )) + 
	geom_line ( data = pr_nspp_log, aes(x = (nspp), y = (Biomass),color = "1" ) ) + 
	geom_point (data= rDIT_lake_eq,aes(x = (shannon), y = (Biomass),color = "2")) +
	geom_line ( data = pr_H, aes(x = shannon, y = Biomass,color = "2" ) ) + 
	geom_point (data= rDIT_lake_eq,aes(x = (rS), y =(Biomass),color = "3")) +
	geom_line ( data = pr_rS, aes(x = rS, y = Biomass,color = "3" ) ) + 
	geom_point( data= rDIT_lake_eq,aes (x = (rMI), y=(Biomass),color = "4" ) ) +
	geom_line ( data = pr_rMI, aes(x = rMI, y = Biomass,color = "4" ) ) + 
	geom_point( data= rDIT_lake_eq,aes (x = (rCE), y=(Biomass),color = "5" ) ) +
	geom_line ( data = pr_rCE, aes(x = rCE, y = Biomass,color = "5" ) ) + 
	scale_y_log10()+ 
	xlab("#Species, Bits")+
	ylab("Biomass")+
	scale_color_discrete(name ="", labels = c("# Species", "SDI", "rSD", "rMI","rCE") )


