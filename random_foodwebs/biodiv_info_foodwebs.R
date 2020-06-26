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
nwebs = 20

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

#Random resources:
c = 0.1
amp = 1
res_R = c(amp,c)

for (w in 1:nwebs){ 
	print(w)

	#Assume 3 trophic levels unless otherwise specified.
	nRsp = ceiling(runif(1)*3)
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

		print( paste( "nRsp", sum(out1[[w]]$out[tl,1:nRsp]>1) ) )
		print( paste( "nCsp", sum(out1[[w]]$out[tl,(nRsp+1):nCsp]>1) ) )
		print( paste( "nPsp", sum(out1[[w]]$out[tl,(nCsp+1):nPsp]>1) ) )		

		# plot(out1[[w]]$out[,1], t="l", ylim = c(0, max(out1[[w]]$out[tl,],na.rm=T) ) )
		# for(n in 2:nRsp){ lines(out1[[w]]$out[,n], col ="red") }
		# for(n in (nRsp+1):(nCsp) ){ lines(out1[[w]]$out[,n], col ="blue") }
		# for(n in (nCsp+1):(nPsp) ){ lines(out1[[w]]$out[,n]) }


	}, error = function(e){}) 


	
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
	
	# rweb1[w] = list(rutledge_web( spp_list=c(nRsp,nCsp,nPsp), pop_ts = out1[[w]]$out[,2:(nspp+1)],
	# 	spp_prms = out1[[w]]$spp_prms) )

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
	di_web[w] = list(get_info_dynamics(pop_ts = floor(out1[[w]]$out[nt1:tl,2:(nspp+1)]), 
		k=k,with_blocks=TRUE))

	## This code takes the population time-series counts output by the ODEs and 
	## calculates the average Transfer Entropy from each species to every other 
	## species. The goal is to get an overview of the major information pathways 
	## in the web.   
	#=============================================================================
	# This function gives:
	# te_web		Average transfer entropy per species as a pairwise matrix
	#=============================================================================
	te_web[w] = list( get_te_web( pop_ts = floor(out1[[w]]$out[nt1:tl,2:(nspp+1)]), 
		k=k) )

	## This code takes the population time-series counts output by the ODEs and 
	## calculates the average Separable Information from each species to every other 
	## species. The goal is to get an overview of the major information pathways 
	## in the web.   
	#=============================================================================
	# This function gives:
	# si_web		Average separable information per species as a pairwise matrix
	#=============================================================================
	si_web[w] = list( get_si_web( pop_ts = floor(out1[[w]]$out[nt1:tl,2:(nspp+1)]), 
		k=k) )
}

save(file = "rand_fwebmod6C.var", out1,  di_web,te_web,si_web)

#=============================================================================
# Examine a particular food web more closely: 
#=============================================================================
library(viridis)
library(fields)
library(igraph)
library(visNetwork)

w=1
#=============================================================================
#Export parameters into csv tables for easier reading. 
#  !!! Make sure to set the name of the excel file below!!!!
#=============================================================================
library(xlsx)

var_load = out1[[w]]$spp_prms[5] #These start at variable 5 and go to 14
write.xlsx(var_load, file="spp_prms_rweb1.xlsx", sheetName="sheet1", row.names=FALSE)
for (n in 6:14){
	var_load = out1[[w]]$spp_prms[n]
	sheet = paste("sheet",n-4, sep='')

	write.xlsx(var_load, file="spp_prms_rweb1.xlsx", sheetName=sheet, append=TRUE,row.names=FALSE)
}

#=============================================================================
#Export the average information theoretic quantities into tables.  
#  !!! Make sure to set the name of the excel file below!!!!
#=============================================================================
library(xlsx)

var_load = di_web[[w]]$ee_means #These start at variable 5 and go to 14
write.xlsx(var_load, file="avg_dit_rweb1.xlsx", sheetName="sheet1", row.names=FALSE)
var_load = di_web[[w]]$ai_means #These start at variable 5 and go to 14
write.xlsx(var_load, file="avg_dit_rweb1.xlsx", sheetName="sheet2",append=TRUE,row.names=FALSE)
var_load = di_web[[w]]$te_means #These start at variable 5 and go to 14
write.xlsx(var_load, file="avg_dit_rweb1.xlsx", sheetName="sheet3",append=TRUE,row.names=FALSE)
var_load = di_web[[w]]$si_means #These start at variable 5 and go to 14
write.xlsx(var_load, file="avg_dit_rweb1.xlsx", sheetName="sheet4",append=TRUE,row.names=FALSE)

#=============================================================================
# Plot each of the average information theoretic metrics as a bar graph
#=============================================================================
fig.name = paste("average_dynamics_rweb1.pdf",sep="")
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


#=============================================================================
# Plot the population dynamics
#=============================================================================
out = out1[[w]]$out
nspp = out1[[w]]$spp_prms$nspp
nRsp = out1[[w]]$spp_prms$nRsp
nCsp = out1[[w]]$spp_prms$nCsp
nPsp = out1[[w]]$spp_prms$nPsp
tl = tend/delta1
par(mfrow=c(3,1))
#Resource species in RED
plot(out[,"1"],t="l",col="red",ylim = c(0,max(out[tl,2:(nRsp+1)],na.rm=T)))
for( n in 1:(nRsp) ) {
lines(out[,paste(n)],t="l",col="red")
}

#Consumer species in BLUE 
plot(out[,paste(nRsp+2)],t="l",col="blue",ylim = c(0,max(out[tl,(nRsp+2):(nRsp+nCsp+1)],na.rm=T)))
for( n in ( (nRsp+1):(nRsp+nCsp) ) ) {
lines(out[,paste(n)],t="l",col="blue")
}

#Predator species in BLACK
plot(out[,paste(nRsp+nCsp+2)],t="l",ylim = c(0,max(out[tl,(nRsp+nCsp+2):(nspp+1)],na.rm=T)))
for( n in ((nRsp+nCsp+1):(nspp) ) ) {
lines(out[3900:4000,paste(n)],t="l")
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
		  		#visSave(file="te_graph1p.html", selfcontained = FALSE, background = "white")
  				visExport( type = "pdf", name = "te_web_biggest_1")


######################################################
# Add information storage (AIS or EE) as a self-loop!#
######################################################
edges_tmp = data.frame(from = c(1:length(spp_use)), to =(1:length(spp_use)),weight =(1:length(spp_use))  )
edges_tmp$value = di_web[[1]]$ai_means[spp_use]
te_visn$edges=rbind(te_visn$edges,edges_tmp)

visNetwork(te_visn$nodes, te_visn$edges) %>%
	visEdges(arrows="to", arrowStrikethrough =FALSE  ) %>%
		visOptions(highlightNearest = list(enabled =TRUE, degree =0) )%>%
		  	visIgraphLayout(layout = "layout_in_circle") %>%
		  		visSave(file="ai_te_graph1.html", selfcontained = FALSE, background = "white")
  				#visExport( type = "pdf", name = "te_web_biggest_1")



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
		  		visSave(file="si_graph1.html", selfcontained = FALSE, background = "white")
  				#visExport( type = "pdf", name = "si_web_biggest_1")


#=============================================================================
# Make combined plots of population and dynamic information metrics with time 
#=============================================================================

#===========================================#
#plot1: Info storage (Excess Entropy or AIS)
#===========================================#

# fig.name = paste("dynamic_info_AIS_rweb1.pdf",sep="")
# pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

#When the figure is only over a subset of the time to show transient dynamics: 
fig.name = paste("dynamic_info_AIS_rweb1_sub.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)


layout.matrix=matrix(c(1:12), nrow = 6, ncol = 2)
layout(mat = layout.matrix,
       heights = c(1.5, 3.5,1.5, 3.5, 1.5, 3.5,
       1.5, 3.5,1.5, 3.5, 1.5, 3.5), # Heights of the rows
       widths = c(12,1)) # Widths of columns

#layout.show(12)

#par(mfrow=c(2,1),mai= c( 0.0, 0.2, 0.0, 0.2), omi=c(0.5,0.75,0.5,0.75)) #,mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))

###Common figure properties

t1 = 5840
nlevel = 64 #For viridis color scheme
#nt_use = dim(di_web[[w]]$ai_local)[1]
nt_use = 5940
rs1 = 450 #lower bound for Resource population plot 

par(oma = c(3,2,3,3) )

#===========================================#
#===========================================#
###Predator species
par( mar = c(0.5,4,0,4) )

plot(out[t1:nt_use,paste(nRsp+nCsp+2)],t="l",ylim = c(0,max(out[t1:nt_use,(nRsp+nCsp+2):(nspp+1)],na.rm=T)), 
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
plot(out[t1:nt_use,paste(nRsp+2)],t="l",col="blue",ylim = c(0,max(out[t1:nt_use,(nRsp+2):(nRsp+nCsp+1)],na.rm=T))
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
plot(out[t1:nt_use,"1"],t="l",col="red",ylim = c(rs1,max(out[t1:nt_use,2:(nRsp+1)],na.rm=T)), ylab="Population", xlab="", 
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

# fig.name = paste("dynamic_info_TE_rweb1.pdf",sep="")
# pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

#When the figure is only over a subset of the time to show transient dynamics: 
fig.name = paste("dynamic_info_TE_rweb1_sub.pdf",sep="")
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
t1 = 5840
nlevel = 64 #For viridis color scheme
#nt_use = dim(di_web[[w]]$ai_local)[1]
nt_use = 5940
rs1 = 450 #lower bound for Resource population plot 

par(oma = c(3,2,3,3) )

#===========================================#
#===========================================#
###Predator species
par( mar = c(0.5,4,0,4) )

plot(out[t1:nt_use,paste(nRsp+nCsp+2)],t="l",ylim = c(0,max(out[t1:nt_use,(nRsp+nCsp+2):(nspp+1)],na.rm=T)), 
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
plot(out[t1:nt_use,paste(nRsp+2)],t="l",col="blue",ylim = c(0,max(out[t1:nt_use,(nRsp+2):(nRsp+nCsp+1)],na.rm=T))
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
plot(out[1:tl,"1"],t="l",col="red",ylim = c(rs1,max(out[t1:nt_use,2:(nRsp+1)],na.rm=T)), ylab="Population", xlab="", 
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
# fig.name = paste("dynamic_info_SI_rweb1.pdf",sep="")
# pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

#When the figure is only over a subset of the time to show transient dynamics: 
fig.name = paste("dynamic_info_SI_rweb1_sub.pdf",sep="")
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
t1 = 5840
nlevel = 64 #For viridis color scheme
#nt_use = dim(di_web[[w]]$ai_local)[1]
nt_use = 5940
rs1 = 450 #lower bound for Resource population plot 

par(oma = c(3,2,3,3) )

#===========================================#
#===========================================#
###Predator species
par( mar = c(0.5,4,0,4) )

plot(out[t1:nt_use,paste(nRsp+nCsp+2)],t="l",ylim = c(0,max(out[tl,(nRsp+nCsp+2):(nspp+1)],na.rm=T)), 
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
plot(out[t1:nt_use,paste(nRsp+2)],t="l",col="blue",ylim = c(0,max(out[t1:nt_use,(nRsp+2):(nRsp+nCsp+1)],na.rm=T))
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
plot(out[t1:nt_use,"1"],t="l",col="red",ylim = c(rs1,max(out[tl,2:(nRsp+1)],na.rm=T)), ylab="Population", xlab="", 
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

