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
# 2. Generate a series of foodwebs building from simple to complex structure
# 3. Perturb the food web in one of two ways: 
#	 A. Remove a species
#	 B. Add a pulse or press perturbation to a species growth rate
# 3. Use information theory to track the resulting food-web structures. 
# 4. This file has a lot of code for visualizing output of both the foodweb 
#	 its information theoretic properties after the main loop. 
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

###Build a series of scenarios going from simple to more complex dynamics
#Number of food webs to generate
nwebs = 6
# scenarios = list(matrix(0,nwebs,1))

###
#Output of each web
out1 = list(matrix(0,nwebs,1))
#Converting the web to Rutledge's compartment model and calculating the information
#theoretic quantities: Shannon Entropy, Mutual Information, Conditional Entropy
rweb1 = list(matrix(0,nwebs,1))
rweb1_eq = rweb1 #Equilibrium dynamics
rweb1_tr = rweb1 #Transient dynamics

#Dynamic information metrics calculated from the (discretized) time series 
di_web = list(matrix(0,nwebs,1))
di_web_eq = di_web #Equilibrium dynamics
di_web_tr = di_web #Transient dynamics

#Track the average transfer entropy and separable information between each pair of 
#species as a way to build a network of information flow through the network. 
te_web = list(matrix(0,nwebs,1))
si_web = list(matrix(0,nwebs,1)) 
te_web_eq = te_web #Equilibrium dynamics
si_web_eq = si_web
te_web_tr = te_web #Transient dynamics
si_web_tr = si_web

#Random resources:
 c = 0
 amp = 1/exp(1)
 res_R = c(amp,c)

#or 
# res_R = NULL

# scenarios[[1]] = list(nRsp = 1, nCsp =0, nPsp = 0)
# scenarios[[2]] = list(nRsp = 3, nCsp =0, nPsp = 0)
# scenarios[[3]] = list(nRsp = 3, nCsp =0, nPsp = 0)
# scenarios[[4]] = list(nRsp = 2, nCsp =1, nPsp = 0)
# scenarios[[5]] = list(nRsp = 1, nCsp =1, nPsp = 1)

#The structure of this code is based on taking an initial food web and 
#going through a series of perturbations. The "w" index now corresponds
#to each of the perturbations. 
w = 1 
#for (w in 1:nwebs){ 
print(w)

#Assume 3 trophic levels unless otherwise specified.
nRsp = 2 #ceiling(runif(1)*30)
nCsp = 2 #ceiling(runif(1)*20)
nPsp = 2 #ceiling(runif(1)*10)
nspp = nRsp+nCsp+nPsp

#Randomly generate the species parameters for the model as well: 
spp_prms = NULL
#Resource: Nearly identical resource dynamics: 
spp_prms$rR = matrix(rnorm(nRsp,1.5,0), nRsp, 1) #intrinsic growth
spp_prms$Ki = matrix(rnorm(nRsp,1.5,0), nRsp, 1) #carrying capacity

#Consumers: 
spp_prms$rC = matrix(0.6, nCsp, 1) #matrix(rnorm(nCsp,.5,0.2), nCsp, 1) #intrisic growth
#spp_prms$rC[2] = 0.199999 #matrix(rnorm(nCsp,.5,0.2), nCsp, 1) #intrisic growth
spp_prms$eFc = matrix(1,nCsp,nRsp) # just make the efficiency for everything 1 for now
spp_prms$muC = matrix(0.3, nCsp, 1) #matrix(rnorm(nCsp,0.6,0.1), nCsp, 1) #mortality rates
#Consumption rates: 
#Generate a hierarchy where each species predominantly feeds on particular resource. 
dspp = abs((nCsp - nRsp))
#hier1= seq(1/nRsp, (1-1/nRsp), length=nRsp)
hier1 = c(matrix(0.9,nRsp,1)) #1
#hier1 = c(matrix(0.2,nRsp,1)) #2

spp_prms$cC = hier1 
for( n in 1:nCsp) {
	spp_prms$cC = cbind(spp_prms$cC, shifter(hier1,n))
}
spp_prms$cC = matrix(spp_prms$cC[1:nRsp,1:nCsp ],nRsp,nCsp)
#spp_prms$cC[,1] = c(0.5,0.1)
#spp_prms$cC[,2] = c(0.5,0.1)

#  spp_prms$cC[,1] = c(0.5,0.0,0.1)
#  spp_prms$cC[,2] = c(0.1,0.6,0.1)
#  spp_prms$cC[,3] = c(0.1,0.0,0.5)

# #Predators: 
spp_prms$rP =  matrix(0.4, nPsp, 1) #matrix(rnorm(nPsp,0.5,0), nPsp, 1) #intrisic growth
spp_prms$eFp = matrix(1,nPsp,nCsp) # just make the efficiency for everything 1 for now
spp_prms$muP = matrix(0.2, nPsp, 1)#matrix(rnorm(nPsp,0.6,0), nPsp, 1)  #mortality rates
#Consumption rates: 
#Generate a hierarchy where each species predominantly feeds on particular resource. 
dspp = ((nPsp - nCsp))
if(dspp<0){dspp = 0 }
#hier1= seq(1/nCsp, (1-1/nCsp), length = nCsp)
hier1 = c(matrix(0.1,nCsp,1)) #1
#hier1 = c(matrix(0.5,nCsp,1)) #2
spp_prms$cP = hier1
for( n in 1:nPsp) {
	spp_prms$cP = cbind(spp_prms$cP, shifter(hier1,n))
}
spp_prms$cP = matrix(spp_prms$cP[1:nCsp,1:nPsp],nCsp,nPsp)
# spp_prms$cP[,1] = c(0.5,0.4)
# spp_prms$cP[,2] = c(0.5,0.4)

 # spp_prms$cP[,1] = c(0.5,0.1,0.1)
 # spp_prms$cP[,2] = c(0.1,0.6,0.1)
 # spp_prms$cP[,3] = c(0.1,0.1,0.5)

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
#winit = runif(nspp,min=1,max=6)
#winit = c(1.007368, 1.007368, 1.005849, 1.005849, 0.9988030, 0.9988030)
tryCatch( {out1[w] = list(food_web_dynamics (spp_list = c(nRsp,nCsp,nPsp), spp_prms = spp_prms, 
	tend, delta1, res_R = res_R,final = FALSE ))}, error = function(e){}) 

out1[[w]]$out[tl,]


#=============================================================================
# Perturbation set 1: Remove each species and track the dynamics
#=============================================================================
for (s in 1:nspp){

	out_temp =NULL
	out_temp2 =NULL
	inv_spp = s
	winit =  out1[[1]]$out[tl,2:(nspp+1)]
	winit[inv_spp] = 0

	#Equilibrate new community
	tryCatch( {out_temp = (food_web_dynamics (spp_list = c(nRsp,nCsp,nPsp), spp_prms = spp_prms, 
	tend, delta1, winit = winit, res_R = res_R,final = FALSE ))}, error = function(e){}) 

	#Invade
	#ti = which(out_temp$out[ ,nRsp+3] == max(out_temp$out[,nRsp+3] ) )
	ti = tl
	winit =  out_temp$out[ti,2:(nspp+1)]
	winit[inv_spp] = .1

	#Invade new community
	tryCatch( {out_temp2 = (food_web_dynamics (spp_list = c(nRsp,nCsp,nPsp), spp_prms = spp_prms, 
	tend, delta1, winit = winit, res_R = res_R,final = FALSE ))}, error = function(e){}) 

	#Competition from competitor
	#plot(out_temp2$out[1500:2000,(s+1)])
	#ts1=log(out_temp2$out[1500:2000,(s+1)])
	
	#"Competition" from predator
	plot(out_temp2$out[50:300,(s+1)])
	ts1=log(out_temp2$out[50:300,(s+1)])
	ts2=1:length(ts1)
	summary(lm(ts1~ts2))


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
	f1 = 100 #scaling factor
	di_web[w] = list(get_info_dynamics(pop_ts = floor(f1*out1[[w]]$out[nt1:tl,2:(nspp+1)]), 
		k=k,with_blocks=FALSE))

	## This code takes the population time-series counts output by the ODEs and 
	## calculates the average Transfer Entropy from each species to every other 
	## species. The goal is to get an overview of the major information pathways 
	## in the web.   
	#=============================================================================
	# This function gives:
	# te_web		Average transfer entropy per species as a pairwise matrix
	#=============================================================================
	te_web[w] = list( get_te_web( pop_ts = floor(f1*out1[[w]]$out[nt1:tl,2:(nspp+1)]), 
		k=k) )

	## This code takes the population time-series counts output by the ODEs and 
	## calculates the average Separable Information from each species to every other 
	## species. The goal is to get an overview of the major information pathways 
	## in the web.   
	#=============================================================================
	# This function gives:
	# si_web		Average separable information per species as a pairwise matrix
	#=============================================================================
	si_web[w] = list( get_si_web( pop_ts = floor(f1*out1[[w]]$out[nt1:tl,2:(nspp+1)]), 
		k=k) )


	#=============================================================================
	# This section performs the same information theoretic calculations as above, 
	# but in a region of dynamics corresponding to (ideally) equilibrium conditions. 
	#=============================================================================
	nt1e = tl/2
	nt2e = tl
	rweb1_eq[w] = list(rutledge_web( spp_list=c(nRsp,nCsp,nPsp), pop_ts = out1[[w]]$out[nt1e:tl,2:(nspp+1)],
	spp_prms = out1[[w]]$spp_prms) )
	di_web_eq[w] = list(get_info_dynamics(pop_ts = floor(out1[[w]]$out[nt1e:tl,2:(nspp+1)]), 
		k=k,with_blocks=TRUE))
	te_web_eq[w] = list( get_te_web( pop_ts = floor(out1[[w]]$out[nt1e:tl,2:(nspp+1)]), 
		k=k) )
	si_web_eq[w] = list( get_si_web( pop_ts = floor(out1[[w]]$out[nt1e:tl,2:(nspp+1)]), 
		k=k) )

	#=============================================================================
	# This section performs the same information theoretic calculations as above, 
	# but in a region of dynamics corresponding to (ideally) transient conditions. 
	#=============================================================================
	nt1t = 1
	nt2t = 1000
	rweb1_tr[w] = list(rutledge_web( spp_list=c(nRsp,nCsp,nPsp), pop_ts = out1[[w]]$out[nt1t:nt2t,2:(nspp+1)],
	spp_prms = out1[[w]]$spp_prms) )
	di_web_tr[w] = list(get_info_dynamics(pop_ts = floor(out1[[w]]$out[nt1t:nt2t,2:(nspp+1)]), 
		k=k,with_blocks=TRUE))
	te_web_tr[w] = list( get_te_web( pop_ts = floor(out1[[w]]$out[nt1t:nt2t,2:(nspp+1)]), 
		k=k) )
	si_web_tr[w] = list( get_si_web( pop_ts = floor(out1[[w]]$out[nt1t:nt2t,2:(nspp+1)]), 
		k=k) )

}

save(file = "scen_fwebmod6Rand.var", out1, rweb1, di_web,te_web,si_web, 
	rweb1_eq, di_web_eq,te_web_eq,si_web_eq, 
	rweb1_tr, di_web_tr,te_web_tr,si_web_tr)

#=============================================================================
# Examine a particular food web more closely: 
#=============================================================================
library(viridis)
library(fields)
library(igraph)
library(visNetwork)

w=1

#=============================================================================
# Plot the population dynamics
#=============================================================================
out = out1[[w]]$out
nspp = out1[[w]]$spp_prms$nspp
nRsp = out1[[w]]$spp_prms$nRsp
nCsp = out1[[w]]$spp_prms$nCsp
nPsp = out1[[w]]$spp_prms$nPsp
#tlg = tend/delta1
tlg = tl

par(mfrow=c(3,1))
#Resource species in RED
plot(out[1:tlg,2],t="l",col="red",ylim = c(0,max(out[,2:(nRsp+1)],na.rm=T)))
for( n in 2:(nRsp+1) ) {
lines(out[1:tlg,n],t="l",col="red")
}

#Consumer species in BLUE 
plot(out[1:tlg,nRsp+2],t="l",col="blue",ylim = c(0,max(out[,(nRsp+2):(nRsp+nCsp+1)],na.rm=T)))
for( n in (nRsp+2):(nRsp+nCsp+1)  ) {
lines(out[1:tlg,n],t="l",col="blue")
}

#Predator species in BLACK
plot(out[1:tlg,paste(nRsp+nCsp+1)],t="l",ylim = c(0,max(out[,(nRsp+nCsp+2):(nspp+1)],na.rm=T)))
for( n in ((nRsp+nCsp+1):(nspp) ) ) {
lines(out[1:tlg,paste(n)],t="l")
}
#=============================================================================
# Network plot of the food web
#=============================================================================
for(w in 1:nwebs) {

fig.name = paste("food_web_test5",w,"_R",nRsp,"_C",nCsp,"_P",nPsp,".html", sep="")
fig.name2 = paste("food_web_test5",w,"_R",nRsp,"_C",nCsp,"_P",nPsp,".pdf", sep="")

nRsp = out1[[w]]$spp_prms$nRsp
nCsp = out1[[w]]$spp_prms$nCsp
nPsp = out1[[w]]$spp_prms$nPsp

nspp = nRsp+nCsp+nPsp
###This shows the network, but only highlights the largest link between each
###node
#Pair down the graph by removing species that have essentially gone extinct
#from the system. 
pop_web1 = matrix(0,nspp,nspp)
pop_web1 [1:nRsp,1:nRsp] = diag(1,nRsp,nRsp)
#pop_web1 [(nRsp+1):(nRsp+nCsp),(nRsp+1):(nRsp+nCsp)] = out1[[w]]$spp_prms$cC
pop_web1 [(nRsp+1):(nRsp+nCsp),(nRsp+1):(nRsp+nCsp)] = 
	diag(c(out1[[w]]$spp_prms$muC),nCsp,nCsp)
pop_web1 [1:nRsp,(nRsp+1):(nRsp+nCsp)] = out1[[w]]$spp_prms$cC

pop_web1 [(nRsp+1):(nRsp+nCsp),(nRsp+nCsp+1):(nspp)] = out1[[w]]$spp_prms$cP
pop_web1 [(nRsp+nCsp+1):(nspp),(nRsp+nCsp+1):(nspp)] = 
	diag(c(out1[[w]]$spp_prms$muP),nPsp,nPsp)

#Make an igraph object
pop_gr = graph_from_adjacency_matrix(pop_web1, mode="directed", weighted=T)
#Convert to VisNetwork list
pop_visn = toVisNetworkData(pop_gr)
pop_visn$nodes$value = pop_visn$nodes$id
#Copy column "weight" to new column "value" in list "edges"
#pop_visn$edges$value = pop_visn$edges$weight
#Color code the nodes by trophic level 
spp_colors= c( matrix("red",nRsp,1),matrix("blue",nCsp,1),
	matrix("black",nPsp,1) )
pop_visn$nodes$color = spp_colors

#Plot this as an HTML object 
#Add arrows to show direction
#Add an option that when a node is clicked on only the "from" arrows are shown
visNetwork(pop_visn$nodes, pop_visn$edges) %>%
	visEdges(arrows="to", arrowStrikethrough =FALSE  ) %>%
		visOptions(highlightNearest = list(enabled =TRUE, degree =0) )%>%
		  	visIgraphLayout(layout = "layout_in_circle") %>%
		  		#visSave(file=fig.name, selfcontained = FALSE, background = "white")
  				visExport( type = "pdf", name = fig.name2)

}
#=============================================================================
#Export parameters into csv tables for easier reading. 
#  !!! Make sure to set the name of the excel file below!!!!
#=============================================================================
library(xlsx)

for(w in 1:nwebs) {
file.name = paste("spp_prms_sweb",w,".xlsx",sep="")

var_load = out1[[w]]$spp_prms[5] #These start at variable 5 and go to 14
write.xlsx(var_load, file=file.name, sheetName="sheet1", row.names=FALSE)
for (n in 6:14){
	var_load = out1[[w]]$spp_prms[n]
	sheet = paste("sheet",n-4, sep='')

	write.xlsx(var_load, file=file.name, sheetName=sheet, append=TRUE,row.names=FALSE)
}

}
#=============================================================================
#Export the average information theoretic quantities into tables.  
#  !!! Make sure to set the name of the excel file below!!!!
#=============================================================================
library(xlsx)
file.name = paste("avg_dit_sweb",w,".xlsx",sep="")

var_load = di_web[[w]]$ee_means #These start at variable 5 and go to 14
write.xlsx(var_load, file=file.name, sheetName="sheet1", row.names=FALSE)
var_load = di_web[[w]]$ai_means #These start at variable 5 and go to 14
write.xlsx(var_load, file=file.name, sheetName="sheet2",append=TRUE,row.names=FALSE)
var_load = di_web[[w]]$te_means #These start at variable 5 and go to 14
write.xlsx(var_load, file=file.name, sheetName="sheet3",append=TRUE,row.names=FALSE)
var_load = di_web[[w]]$si_means #These start at variable 5 and go to 14
write.xlsx(var_load, file=file.name, sheetName="sheet4",append=TRUE,row.names=FALSE)

#=============================================================================
# Plot each of the average information theoretic metrics as a bar graph
#=============================================================================

for(w in 1:nwebs) {

fig.name = paste("average_dynamics_sweb_test6",w,".pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

layout.matrix=matrix(c(1:4), nrow = 2, ncol = 2)
layout(mat = layout.matrix,
       heights = c(5,5), # Heights of the rows
       widths = c(5,5)) # Widths of columns

#layout.show(4)

barplot(di_web[[w]]$ee_means,cex.lab =1.3, beside = TRUE,ylab="Bits of information", xlab = "")
abline(v =out1[[w]]$spp_prms$nRsp+1,col="red")
mtext("Resour", side=1, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(v =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp +1,col="blue" )
mtext("Consum", side=1, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Pred", side=1, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )

barplot(di_web[[w]]$ai_means,cex.lab =1.3, beside = TRUE,ylab="Bits of information", xlab = "Species #")
abline(v =out1[[w]]$spp_prms$nRsp+1,col="red" )
mtext("Resour", side=1, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(v =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp +1,col="blue" )
mtext("Consum", side=1, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Pred", side=1, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )
mtext("Average Information Storage", side = 3, line =4)


barplot(di_web[[w]]$te_means,cex.lab =1.3, beside = TRUE,ylab="", xlab = "")
abline(v =out1[[w]]$spp_prms$nRsp+1,col="red" )
mtext("Resour", side=1, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(v =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp +1,col="blue"  )
mtext("Consum", side=1, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Pred", side=1, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )
mtext("Average Information Transfer", side = 3, line = 2)


barplot(di_web[[w]]$si_means,cex.lab =1.3, beside = TRUE,ylab="", xlab = "Species #")
abline(v =out1[[w]]$spp_prms$nRsp+1,col="red"  )
mtext("Resour", side=1, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(v =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp +1 ,col="blue" )
mtext("Consum", side=1, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Pred", side=1, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )
mtext("Average Information Modification", side = 3, line = 2)

dev.off()

}

#=============================================================================
# Equilibrium region only: 
# Plot each of the average information theoretic metrics as a bar graph
#=============================================================================

for(w in 1:nwebs) {

#fig.name = paste("average_dynamics_sweb_eq",w,".pdf",sep="")
#pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

fig.name = paste("average_dynamics_sweb_eq",w,".png",sep="")
png(file=fig.name, height=5, width=5, units = "in",res=300)



layout.matrix=matrix(c(1:4), nrow = 2, ncol = 2)
layout(mat = layout.matrix,
       heights = c(5,5), # Heights of the rows
       widths = c(5,5)) # Widths of columns

#layout.show(4)

barplot(di_web_eq[[w]]$ee_means,cex.lab =1.3, beside = TRUE,ylab="Bits of information", xlab = "")
abline(v =out1[[w]]$spp_prms$nRsp+1,col="red")
mtext("Resour", side=1, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(v =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp +1,col="blue" )
mtext("Consum", side=1, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Pred", side=1, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )

barplot(di_web_eq[[w]]$ai_means,cex.lab =1.3, beside = TRUE,ylab="Bits of information", xlab = "Species #")
abline(v =out1[[w]]$spp_prms$nRsp+1,col="red" )
mtext("Resour", side=1, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(v =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp +1,col="blue" )
mtext("Consum", side=1, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Pred", side=1, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )
mtext("Average Information Storage", side = 3, line =4)


barplot(di_web_eq[[w]]$te_means,cex.lab =1.3, beside = TRUE,ylab="", xlab = "")
abline(v =out1[[w]]$spp_prms$nRsp+1,col="red" )
mtext("Resour", side=1, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(v =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp +1,col="blue"  )
mtext("Consum", side=1, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Pred", side=1, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )
mtext("Average Information Transfer", side = 3, line = 2)


barplot(di_web_eq[[w]]$si_means,cex.lab =1.3, beside = TRUE,ylab="", xlab = "Species #")
abline(v =out1[[w]]$spp_prms$nRsp+1,col="red"  )
mtext("Resour", side=1, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(v =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp +1 ,col="blue" )
mtext("Consum", side=1, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Pred", side=1, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )
mtext("Average Information Modification", side = 3, line = 2)

dev.off()

}

#=============================================================================
# Transient dynamics region only: 
# Plot each of the average information theoretic metrics as a bar graph
#=============================================================================
for(w in 1:nwebs) {

#fig.name = paste("average_dynamics_sweb_tr",w,".pdf",sep="")
#pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

fig.name = paste("average_dynamics_sweb_tr",w,".png",sep="")
png(file=fig.name, height=5, width=5, units = "in",res=300)


layout.matrix=matrix(c(1:4), nrow = 2, ncol = 2)
layout(mat = layout.matrix,
       heights = c(5,5), # Heights of the rows
       widths = c(5,5)) # Widths of columns

#layout.show(4)

barplot(di_web_tr[[w]]$ee_means,cex.lab =1.3, beside = TRUE,ylab="Bits of information", xlab = "")
abline(v =out1[[w]]$spp_prms$nRsp+1,col="red")
mtext("Resour", side=1, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(v =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp +1,col="blue" )
mtext("Consum", side=1, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Pred", side=1, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )

barplot(di_web_tr[[w]]$ai_means,cex.lab =1.3, beside = TRUE,ylab="Bits of information", xlab = "Species #")
abline(v =out1[[w]]$spp_prms$nRsp+1,col="red" )
mtext("Resour", side=1, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(v =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp +1,col="blue" )
mtext("Consum", side=1, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Pred", side=1, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )
mtext("Average Information Storage", side = 3, line =4)


barplot(di_web_tr[[w]]$te_means,cex.lab =1.3, beside = TRUE,ylab="", xlab = "")
abline(v =out1[[w]]$spp_prms$nRsp+1,col="red" )
mtext("Resour", side=1, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(v =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp +1,col="blue"  )
mtext("Consum", side=1, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Pred", side=1, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )
mtext("Average Information Transfer", side = 3, line = 2)


barplot(di_web_tr[[w]]$si_means,cex.lab =1.3, beside = TRUE,ylab="", xlab = "Species #")
abline(v =out1[[w]]$spp_prms$nRsp+1,col="red"  )
mtext("Resour", side=1, at = c( out1[[w]]$spp_prms$nRsp/2 ) )
abline(v =out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp +1 ,col="blue" )
mtext("Consum", side=1, at = c( (out1[[w]]$spp_prms$nCsp+out1[[w]]$spp_prms$nRsp )-(out1[[w]]$spp_prms$nCsp)/2 ) )
mtext("Pred", side=1, at = c( nspp-(out1[[w]]$spp_prms$nPsp)/2 ) )
mtext("Average Information Modification", side = 3, line = 2)

dev.off()

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


#=============================================================================
# Network plots of information storage.
#	The local exess entropy or active information could be used to show the 
#	dominant cycles involved in information storage...
#=============================================================================


#=============================================================================
# Network plots of information transfer.
# 	This uses the average Transfer Entropy between each species pair to create
#	a directed network of information transfers. 
#=============================================================================

for(w in 1:nwebs) {

#fig.name = paste("average_dynamics_sweb_eq",w,".png",sep="")
#png(file=fig.name, height=5, width=5, units = "in",res=300)

fig.name = paste("te_graph_test6",w,"_R",nRsp,"_C",nCsp,"_P",nPsp,".html", sep="")

###This shows the network, but only highlights the largest link between each
###node
#Pair down the graph by removing species that have essentially gone extinct
#from the system. 
spp_use = (1:nspp)[out1[[w]]$out[10000,2:nspp]>1e-5]
te_web1 = te_web[[w]][spp_use,spp_use]
#Make an igraph object
te_gr = graph_from_adjacency_matrix(te_web1, mode="directed", weighted=T)
#Convert to VisNetwork list
te_visn = toVisNetworkData(te_gr)
te_visn$nodes$value = te_visn$nodes$id
#Copy column "weight" to new column "value" in list "edges"
te_visn$edges$value = te_visn$edges$weight
#Further prune links that are smaller than the 95% interval
m1 = mean(c(log(te_visn$edges$value)))
sd1 = sqrt(var(c(log(te_visn$edges$value))))
te_visn$edges =te_visn$edges[log(te_visn$edges$value) > (m1-sd1), ]

#Color code the nodes by trophic level 
spp_colors= c( matrix("red",nRsp,1),matrix("blue",nCsp,1),
	matrix("black",nPsp,1) )
spp_colors = spp_colors [spp_use]
te_visn$nodes$color = spp_colors


#Plot this as an HTML object 
#Add arrows to show direction
#Add an option that when a node is clicked on only the "from" arrows are shown
visNetwork(te_visn$nodes, te_visn$edges) %>%
	visEdges(arrows="to", arrowStrikethrough =FALSE  ) %>%
		visOptions(highlightNearest = list(enabled =TRUE, degree =0) )%>%
		  	visIgraphLayout(layout = "layout_in_circle") %>%
		  		visSave(file=fig.name, selfcontained = FALSE, background = "white")
  				#visExport( type = "pdf", name = fig.name)
}

######################################################
# Add information storage (AIS or EE) as a self-loop!#
######################################################
for(w in 1:nwebs) {

fig.name = paste("ai_te_graph_test6",w,"_R",nRsp,"_C",nCsp,"_P",nPsp,".html", sep="")

edges_tmp = data.frame(from = c(1:length(spp_use)), to =(1:length(spp_use)),weight =(1:length(spp_use))  )
edges_tmp$value = di_web[[1]]$ai_means[spp_use]
te_visn$edges=rbind(te_visn$edges,edges_tmp)

visNetwork(te_visn$nodes, te_visn$edges) %>%
	visEdges(arrows="to", arrowStrikethrough =FALSE  ) %>%
		visOptions(highlightNearest = list(enabled =TRUE, degree =0) )%>%
		  	visIgraphLayout(layout = "layout_in_circle") %>%
		  		visSave(file=fig.name, selfcontained = FALSE, background = "white")
  				#visExport( type = "pdf", name = "te_web_biggest_1")

}

######################################
#Because transfer can be asymmetrical, make 2 different graphs showing direction
#of flows. 
te_gr1 = graph_from_adjacency_matrix( (te_web[[w]]*lower.tri(te_web[[w]])), mode="directed", weighted=T)
te_gr2 = graph_from_adjacency_matrix( (te_web[[w]]*upper.tri(te_web[[w]])), mode="directed", weighted=T)

#Convert to VisNetwork list
te_visn1 = toVisNetworkData(te_gr1)
te_visn2 = toVisNetworkData(te_gr2)
#Copy column "weight" to new column "value" in list "edges"
te_visn1$edges$value = te_visn1$edges$weight
te_visn2$edges$value = te_visn2$edges$weight
#Color code the nodes by trophic level 
te_visn1$nodes$color = c( matrix("red",nRsp,1),matrix("blue",nCsp,1),
	matrix("black",nPsp,1) )
te_visn2$nodes$color = c( matrix("red",nRsp,1),matrix("blue",nCsp,1),
	matrix("black",nPsp,1) )
#te_visn1$nodes$color = c( matrix(c("red","blue","black"),9,1) )
#te_visn2$nodes$color = c( matrix(c("red","blue","black"),9,1) )

#Plot this as an HTML object 
#Add arrows to show direction: 
te_visn1$edges$arrows = c(matrix("to",dim(te_visn1$edges)[1]))
te_visn2$edges$arrows = c(matrix("to",dim(te_visn2$edges)[1]))

visNetwork(te_visn1$nodes, te_visn1$edges) %>%
  visIgraphLayout(layout = "layout_in_circle") %>%
  	visExport( type = "pdf", name = "te_web_clock_1")

visNetwork(te_visn2$nodes, te_visn2$edges) %>%
  visIgraphLayout(layout = "layout_in_circle") %>%
  	visExport( type = "pdf", name = "te_web_clock_1")


#=============================================================================
# Network plots of information modification.
# 	This uses the average Separable Information between each species pair to create
#	a directed network of information transfers. 
#=============================================================================
for(w in 1:nwebs) {

fig.name = paste("si_graph",w,"_R",nRsp,"_C",nCsp,"_P",nPsp,".html", sep="")

###This shows the network, but only highlights the largest link between each
###node
#Pair down the graph by removing species that have essentially gone extinct
#from the system. 
spp_use = (1:nspp)[out1[[w]]$out[10000,2:nspp]>1e-5]
si_web1 = si_web[[w]][spp_use,spp_use]
#Make an igraph object
si_gr = graph_from_adjacency_matrix(si_web1, mode="directed", weighted=T)
#Convert to VisNetwork list
si_visn = toVisNetworkData(si_gr)
si_visn$nodes$value = si_visn$nodes$id
#Copy column "weight" to new column "value" in list "edges"
si_visn$edges$value = si_visn$edges$weight
#Color code the nodes by trophic level 
spp_colors= c( matrix("red",nRsp,1),matrix("blue",nCsp,1),
	matrix("black",nPsp,1) )
spp_colors = spp_colors [spp_use]
si_visn$nodes$color = spp_colors

#Plot this as an HTML object 
#Add arrows to show direction
#Add an option that when a node is clicked on only the "from" arrows are shown
visNetwork(si_visn$nodes, si_visn$edges) %>%
	visEdges(arrows="to", arrowStrikethrough =FALSE  ) %>%
		visOptions(highlightNearest = list(enabled =TRUE, degree =0) )%>%
		  	visIgraphLayout(layout = "layout_in_circle") %>%
		  		visSave(file=fig.name, selfcontained = FALSE, background = "white")
  				#visExport( type = "pdf", name = "si_web_biggest_1")

}
#=============================================================================
# Make combined plots of population and dynamic information metrics with time 
#=============================================================================
#=============================================================================
# Plot the population dynamics
#=============================================================================
#Make a combined eq/invasion plot from other runs if needed: 
#
#

tiA=19000
tiB=1
teA=19990
teB=500

outA = out1[[2]]$out[tiA:teA,]
outB = out1[[3]]$out[tiB:teB,]

diA = di_web[[2]]$ai_local[tiA:teA,]
diB = di_web[[3]]$ai_local[tiB:teB,]
 
trA = di_web[[2]]$te_local[tiA:teA,]
trB = di_web[[3]]$te_local[tiB:teB,]

out1[[7]] = NULL
di_web[[7]] = NULL
out1[[7]]= list(out=rbind(outA,outB),spp_prms=NULL)
di_web[[7]] = list(ai_local=rbind(diA,diB), te_local = rbind(trA,trB) )
#=============================================================================


out = out1[[w]]$out
nspp = out1[[w]]$spp_prms$nspp
nRsp = out1[[w]]$spp_prms$nRsp
nCsp = out1[[w]]$spp_prms$nCsp
nPsp = out1[[w]]$spp_prms$nPsp
#tlg = tend/delta1
#tlg = tl
tlg = dim(out)[1]

par(mfrow=c(3,1))
#Resource species in RED
plot(out[1:tlg,2],t="l",col="red",ylim = c(0,max(out[,2:(nRsp+1)],na.rm=T)))
for( n in 2:(nRsp+1) ) {
lines(out[1:tlg,n],t="l",col="red")
}

#Consumer species in BLUE 
plot(out[1:tlg,nRsp+2],t="l",col="blue",ylim = c(0,max(out[,(nRsp+2):(nRsp+nCsp+1)],na.rm=T)))
for( n in (nRsp+2):(nRsp+nCsp+1)  ) {
lines(out[1:tlg,n],t="l",col="blue")
}

#Predator species in BLACK
plot(out[1:tlg,paste(nRsp+nCsp+1)],t="l",ylim = c(0,max(out[,(nRsp+nCsp+2):(nspp+1)],na.rm=T)))
for( n in ((nRsp+nCsp+1):(nspp) ) ) {
lines(out[1:tlg,paste(n)],t="l")
}

#===========================================#
#plot1: Info storage (Excess Entropy or AIS)
#===========================================#

# fig.name = paste("dynamic_info_AIS_sweb1.pdf",sep="")
# pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

#When the figure is only over a subset of the time to show transient dynamics: 
fig.name = paste("dynamic_info_AIS_sweb1test64_sub4.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)


layout.matrix=matrix(c(1:12), nrow = 6, ncol = 2)
layout(mat = layout.matrix,
       heights = c(1.5, 3.5,1.5, 3.5, 1.5, 3.5,
       1.5, 3.5,1.5, 3.5, 1.5, 3.5), # Heights of the rows
       widths = c(12,1)) # Widths of columns

#layout.show(12)

#par(mfrow=c(2,1),mai= c( 0.0, 0.2, 0.0, 0.2), omi=c(0.5,0.75,0.5,0.75)) #,mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))

###Common figure properties

t1=1
#t1 = 14100
nlevel = 64 #For viridis color scheme
#nt_use = dim(di_web[[w]]$ai_local)[1]
#nt_use = tl -100
nt_use = 10000
rs1 = 450 #lower bound for Resource population plot 

par(oma = c(3,2,3,3) )

#===========================================#
#===========================================#
###Predator species
par( mar = c(0.5,4,0,4) )

plot(out[t1:nt_use,paste(nRsp+nCsp+2)],t="l",ylim = c(0,max(out[t1:nt_use,(nRsp+nCsp+2):(nspp+1)],na.rm=T)+1), 
	ylab="Population", xlab="", xaxs="i", xaxt="n",yaxs="i",cex.main=1.2,cex.lab=1.2)
for( n in ((nRsp+nCsp+1):(nspp) ) ) {
lines(out[t1:nt_use,paste(n)],t="l")
}
mtext("Local Information Storage", side = 3, line = 0, outer = TRUE)

#Local excess entropy
#par( mar = c(2,4,0,4) )
# nt_use = dim(di_web[[w]]$ee_local)[1]
#image( 1:nt_use, 1:nCsp, di_web[[w]]$ee_local[,(nRsp+nCsp+1):(nspp)], ylab="Species number", 
#	xlab="Time",col=viridis(nlevel) )

#Local active information storage
par( mar = c(2,4,0,4) )
image( t1:nt_use, 1:nPsp, di_web[[w]]$ai_local[t1:nt_use,(nRsp+nCsp+1):(nspp)], ylab="Species #", 
	xlab="Time",col=viridis(nlevel),cex.main=1.3,cex.lab=1.3)

###Consumer species
par( mar = c(0.5,4,0,4) )
#Consumer species in BLUE
plot(out[t1:nt_use,paste(nRsp+2)],t="l",col="blue",ylim = c(0,max(out[t1:nt_use,(nRsp+2):(nRsp+nCsp+1)],na.rm=T)+1)
	, ylab="Population", xlab="", xaxs="i", xaxt="n",yaxs="i",cex.main=1.2,cex.lab=1.2)
for( n in ( (nRsp+1):(nRsp+nCsp) ) ) {
lines(out[t1:nt_use,paste(n)],t="l",col="blue")
}

#Local excess entropy
#par( mar = c(2,4,0,4) )
# nt_use = dim(di_web[[w]]$ee_local)[1]
#image( 1:nt_use, 1:nCsp, di_web[[w]]$ee_local[,(nRsp+1):(nRsp+nCsp)], ylab="Species number", 
#	xlab="Time",col=viridis(nlevel) )

#Local active information storage
par( mar = c(2,4,0,4) )
image( t1:nt_use, 1:nCsp, di_web[[w]]$ai_local[t1:nt_use,(nRsp+1):(nRsp+nCsp)], ylab="Species #", 
	xlab="Time",col=viridis(nlevel),,cex.main=1.3,cex.lab=1.3 )

###Resource Species
par( mar = c(0.5,4,0,4) )
#Resource species in RED
plot(out[t1:nt_use,"1"],t="l",col="red",ylim = c(0,max(out[t1:nt_use,2:(nRsp+1)],na.rm=T)+1), ylab="Population", xlab="", 
  xaxs="i", xaxt="n",yaxs="i",cex.main=1.2,cex.lab=1.2, )
for( n in 1:(nRsp) ) {
lines(out[t1:nt_use,paste(n)],t="l",col="red")
}

#Local excess entropy
#par( mar = c(2,4,0,4) )
# nt_use = dim(di_web[[w]]$ee_local)[1]
#image( 1:nt_use, 1:nRsp, di_web[[w]]$ee_local[,1:nRsp], ylab="Species number", 
#	xlab="Time",col=viridis(nlevel) )

#Local active information storage
par( mar = c(2,4,0,4) )
image( t1:nt_use, 1:nRsp, di_web[[w]]$ai_local[t1:nt_use,1:nRsp], ylab="Species #", 
	xlab="Time",col=viridis(nlevel),cex.main=1.3,cex.lab=1.3 )

###Plot color bars for image plots: 

#Color bar 1
par( mar = c(0.5,0.5,0.5,0.5) )
frame()
par( mar = c(3,0,0,2) )
var_dist =  di_web[[w]]$ai_local[t1:nt_use,(nRsp+nCsp+1):(nspp)]
image(1,(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	t(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	ylab="",xaxt='n',col=viridis(nlevel))

#Color bar 2
par( mar = c(0.5,0.5,0.5,0.5) )
frame()
par( mar = c(3,0,0,2) )
var_dist =  di_web[[w]]$ai_local[t1:nt_use,(nRsp+nCsp+1):(nspp)]
image(1,(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	t(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	ylab="",xaxt='n',col=viridis(nlevel))


#Color bar 3
par( mar = c(0.5,0.5,0.5,0.5) )
frame()
par( mar = c(3,0,0,2) )
var_dist =  di_web[[w]]$ai_local[t1:nt_use,(nRsp+nCsp+1):(nspp)]
image(1,(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	t(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	ylab="",xaxt='n',col=viridis(nlevel))



dev.off()

#===========================================#
#plot2: Information transmission (TE)
#===========================================#

# fig.name = paste("dynamic_info_TE_sweb1.pdf",sep="")
# pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

#When the figure is only over a subset of the time to show transient dynamics: 
fig.name = paste("dynamic_info_TE_sweb1test64_sub1.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)


layout.matrix=matrix(c(1:12), nrow = 6, ncol = 2)
layout(mat = layout.matrix,
       heights = c(1.5, 3.5,1.5, 3.5, 1.5, 3.5,
       1.5, 3.5,1.5, 3.5, 1.5, 3.5), # Heights of the rows
       widths = c(12,1)) # Widths of columns

#layout.show(12)

#par(mfrow=c(2,1),mai= c( 0.0, 0.2, 0.0, 0.2), omi=c(0.5,0.75,0.5,0.75)) #,mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))

###Common figure properties
nlevel = 64 #For viridis color 
t1=1
#t1 = 14100
nlevel = 64 #For viridis color scheme
#nt_use = dim(di_web[[w]]$ai_local)[1]
nt_use = tl-100
#nt_use = 10000
rs1 = 450 #lower bound for Resource population plot 

par(oma = c(3,2,3,3) )

#===========================================#
#===========================================#
###Predator species
par( mar = c(0.5,4,0,4) )

plot(out[t1:nt_use,paste(nRsp+nCsp+2)],t="l",ylim = c(0,max(out[t1:nt_use,(nRsp+nCsp+2):(nspp+1)],na.rm=T)+1), 
	ylab="Population", xlab="", xaxs="i", xaxt="n",yaxs="i",cex.main=1.2,cex.lab=1.2)
for( n in ((nRsp+nCsp+1):(nspp) ) ) {
lines(out[t1:nt_use,paste(n)],t="l")
}
mtext("Local Transfer Entropy", side = 3, line = 0, outer = TRUE)

#Local Transfer Entropy 
par( mar = c(2,4,0,4) )
image( t1:nt_use, 1:nPsp, di_web[[w]]$te_local[t1:nt_use,(nRsp+nCsp+1):(nspp)], ylab="Species #", 
	xlab="Time",col=viridis(nlevel),cex.main=1.3,cex.lab=1.3)

###Consumer species
par( mar = c(0.5,4,0,4) )
#Consumer species in BLUE
plot(out[t1:nt_use,paste(nRsp+2)],t="l",col="blue",ylim = c(0,max(out[t1:nt_use,(nRsp+2):(nRsp+nCsp+1)],na.rm=T)+1)
	, ylab="Population", xlab="", xaxs="i", xaxt="n",yaxs="i",cex.main=1.2,cex.lab=1.2)
for( n in ( (nRsp+1):(nRsp+nCsp) ) ) {
lines(out[t1:nt_use,paste(n)],t="l",col="blue")
}

#Local transfer entropy
par( mar = c(2,4,0,4) )
image( t1:nt_use, 1:nCsp, di_web[[w]]$te_local[t1:nt_use,(nRsp+1):(nRsp+nCsp)], ylab="Species #", 
	xlab="Time",col=viridis(nlevel),,cex.main=1.3,cex.lab=1.3 )

###Resource Species
par( mar = c(0.5,4,0,4) )
#Resource species in RED

plot(out[t1:nt_use,"1"],t="l",col="red",ylim = c(0,max(out[t1:nt_use,2:(nRsp+1)],na.rm=T)+1), ylab="Population", xlab="", 
  xaxs="i", xaxt="n",yaxs="i",cex.main=1.2,cex.lab=1.2, )
for( n in 1:(nRsp) ) {
lines(out[t1:nt_use,paste(n)],t="l",col="red")
}


#local transfer entropy
par( mar = c(2,4,0,4) )
image( t1:nt_use, 1:nRsp, di_web[[w]]$te_local[t1:nt_use,1:nRsp], ylab="Species #", 
	xlab="Time",col=viridis(nlevel),cex.main=1.3,cex.lab=1.3 )

###Plot color bars for image plots: 

#Color bar 1
par( mar = c(0.5,0.5,0.5,0.5) )
frame()
par( mar = c(3,0,0,2) )
var_dist =  di_web[[w]]$te_local[t1:nt_use,(nRsp+nCsp+1):(nspp)]
image(1,(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	t(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	ylab="",xaxt='n',col=viridis(nlevel))

#Color bar 2
par( mar = c(0.5,0.5,0.5,0.5) )
frame()
par( mar = c(3,0,0,2) )
var_dist =  di_web[[w]]$te_local[t1:nt_use,(nRsp+nCsp+1):(nspp)]
image(1,(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	t(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	ylab="",xaxt='n',col=viridis(nlevel))


#Color bar 3
par( mar = c(0.5,0.5,0.5,0.5) )
frame()
par( mar = c(3,0,0,2) )
var_dist =  di_web[[w]]$te_local[t1:nt_use,(nRsp+nCsp+1):(nspp)]
image(1,(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	t(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	ylab="",xaxt='n',col=viridis(nlevel))

dev.off()

#===========================================#
#plot3: Information modification (SI)
#===========================================#
# fig.name = paste("dynamic_info_SI_sweb1.pdf",sep="")
# pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

#When the figure is only over a subset of the time to show transient dynamics: 
fig.name = paste("dynamic_info_SI_sweb1test54_sub.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)


layout.matrix=matrix(c(1:12), nrow = 6, ncol = 2)
layout(mat = layout.matrix,
       heights = c(1.5, 3.5,1.5, 3.5, 1.5, 3.5,
       1.5, 3.5,1.5, 3.5, 1.5, 3.5), # Heights of the rows
       widths = c(12,1)) # Widths of columns

#layout.show(12)

#par(mfrow=c(2,1),mai= c( 0.0, 0.2, 0.0, 0.2), omi=c(0.5,0.75,0.5,0.75)) #,mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))

###Common figure properties
nlevel = 64 #For viridis color 
t1 = 1
nlevel = 64 #For viridis color scheme
#nt_use = dim(di_web[[w]]$ai_local)[1]
nt_use = tl-100
rs1 = 450 #lower bound for Resource population plot 

par(oma = c(3,2,3,3) )

#===========================================#
#===========================================#
###Predator species
par( mar = c(0.5,4,0,4) )

plot(out[t1:nt_use,paste(nRsp+nCsp+2)],t="l",ylim = c(0,max(out[tl,(nRsp+nCsp+2):(nspp+1)],na.rm=T)+1), 
	ylab="Population", xlab="", xaxs="i", xaxt="n",yaxs="i",cex.main=1.2,cex.lab=1.2)
for( n in ((nRsp+nCsp+1):(nspp) ) ) {
lines(out[t1:nt_use,paste(n)],t="l")
}
mtext("Local Seprable Information", side = 3, line = 0, outer = TRUE)

#Local seprable informatio
par( mar = c(2,4,0,4) )
image( t1:nt_use, 1:nPsp, di_web[[w]]$si_local[t1:nt_use,(nRsp+nCsp+1):(nspp)], ylab="Species #", 
	xlab="Time",col=viridis(nlevel),cex.main=1.3,cex.lab=1.3)

###Consumer species
par( mar = c(0.5,4,0,4) )
#Consumer species in BLUE
plot(out[t1:nt_use,paste(nRsp+2)],t="l",col="blue",ylim = c(0,max(out[t1:nt_use,(nRsp+2):(nRsp+nCsp+1)],na.rm=T)+1)
	, ylab="Population", xlab="", xaxs="i", xaxt="n",yaxs="i",cex.main=1.2,cex.lab=1.2)
for( n in ( (nRsp+1):(nRsp+nCsp) ) ) {
lines(out[t1:nt_use,paste(n)],t="l",col="blue")
}

#Local separable information
par( mar = c(2,4,0,4) )
image( t1:nt_use, 1:nCsp, di_web[[w]]$si_local[t1:nt_use,(nRsp+1):(nRsp+nCsp)], ylab="Species #", 
	xlab="Time",col=viridis(nlevel),,cex.main=1.3,cex.lab=1.3 )

###Resource Species
par( mar = c(0.5,4,0,4) )
#Resource species in RED
plot(out[t1:nt_use,"1"],t="l",col="red",ylim = c(0,max(out[tl,2:(nRsp+1)],na.rm=T)+1), ylab="Population", xlab="", 
  xaxs="i", xaxt="n",yaxs="i",cex.main=1.2,cex.lab=1.2, )
for( n in 1:(nRsp) ) {
lines(out[t1:nt_use,paste(n)],t="l",col="red")
}

#Local separable information
par( mar = c(2,4,0,4) )
image( t1:nt_use, 1:nRsp, di_web[[w]]$si_local[t1:nt_use,1:nRsp], ylab="Species #", 
	xlab="Time",col=viridis(nlevel),cex.main=1.3,cex.lab=1.3 )

###Plot color bars for image plots: 

#Color bar 1
par( mar = c(0.5,0.5,0.5,0.5) )
frame()
par( mar = c(3,0,0,2) )
var_dist =  di_web[[w]]$si_local[t1:nt_use,(nRsp+nCsp+1):(nspp)]
image(1,(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	t(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	ylab="",xaxt='n',col=viridis(nlevel))

#Color bar 2
par( mar = c(0.5,0.5,0.5,0.5) )
frame()
par( mar = c(3,0,0,2) )
var_dist =  di_web[[w]]$si_local[t1:nt_use,(nRsp+nCsp+1):(nspp)]
image(1,(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	t(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	ylab="",xaxt='n',col=viridis(nlevel))


#Color bar 3
par( mar = c(0.5,0.5,0.5,0.5) )
frame()
par( mar = c(3,0,0,2) )
var_dist =  di_web[[w]]$si_local[t1:nt_use,(nRsp+nCsp+1):(nspp)]
image(1,(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	t(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), 
	ylab="",xaxt='n',col=viridis(nlevel))

dev.off()

#===============================================================================
#===============================================================================
#===============================================================================


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

