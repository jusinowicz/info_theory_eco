#=============================================================================
#Load libraries
#=============================================================================
library(tidyverse)
library(RandomFields)
library(vegan)

#Patrick Thompson's metacommunity model functions: 
source("./mcomsimr-master/R/MC_simulate.R")
source("./mcomsimr-master/R/sim_setup.R")
source("../info_theory_functions/info_theory_functions.R")

#=============================================================================
#=============================================================================
# I. Run the metacommunity models! 
#=============================================================================

###This is an example of just running the model. 
# a1=simulate_MC(patches = 16, species = 1, env1Scale = 100, dispersal = 0.001, min_inter = 1, max_inter = 1, env_niche_breadth = 1)
# env.df = a1$env.df
# g<-ggplot(env.df, aes(x = time, y = env1, group = patch, color = factor(patch)))+
#       geom_line()+
#       scale_color_viridis_d(guide=F)
# print(g)

############################################################
### Remove any temporal variance in the landscape 	
### (cannot figure out how to do this with RandomFields): 
############################################################

patches = 16
nspp = 2
landscape = landscape_generate(patches = patches, plot = TRUE)

var_b = 1
env1Scale = 100
timesteps = 1200 
burn_in = 1200  
initialization = 200

#Make the env.df outside of the simulate_MC function with RandomFields.
#This code makes a single time realization of the RandomFields model, then just repeats it over the right number of 
#time steps (timesteps + burn_in)
model = RMexp(var=var_b, scale=env1Scale) + # with variance 4 and scale 10
       RMnugget(var=0) + # nugget
        RMtrend(mean=0.01) # and mean
RF = RFsimulate(model = model, x = landscape$x*10, y = landscape$y*10, spConform=FALSE)
env.df = data.frame(env1 = decostand(RF,method = "range"), patch = 1:nrow(landscape), time = rep(1:(timesteps+burn_in), each = nrow(landscape)))
env.initial = env.df[env.df$time == 1,]

#I've added the "constant = TRUE" option to the function to remove demographic stochasticity. 
#Added mortality to the model. With m = 0, this is the original model. This was necessary to get extinction without
#stochasticity. 
mcm1 =simulate_MC(patches = patches , species = nspp, env1Scale = env1Scale, dispersal = .000001, min_inter = 1.1, max_inter = 1.1,
	 env_niche_breadth = 1, max_r = c(5,4), m = 0.5, env.df = env.df, landscape = landscape, constant = TRUE)

#Plot the environment by patch 
g<-ggplot(mcm1$dynamics.df, aes(x = time, y = env, group = patch, color = factor(patch)))+
      geom_line()+
      scale_color_viridis_d(guide=F)
print(g)
############################################################

#=============================================================================
#=============================================================================
# II. Get the information theoretic metrics! 
# 		There are two levels of metrics: a per-species environmental 
#		information, and a community-level environmental information. 
#=============================================================================
###Environmental information per species: 

#1. Make a time series out of the environmental state and the population. 
#	This is to get the per-species environmental information: 
series1 = data.frame(cbind (env = as.factor(mcm1$dynamics.df$env), species = mcm1$dynamics.df$species, N = mcm1$dynamics.df$N) ,
						time = mcm1$dynamics.df$time) 
series1 = series1[series1$time>=1,]

for (s in 1:nspp){ 

	for (t in 1:(timesteps)) {
		#s1_temp = subset(series1, species == s & time == t)
		s1_temp = spread( subset(series1, time == t), species, N)
		#s1_temp =  subset(series1, time == t)
		s1_temp = s1_temp[,colnames(s1_temp)!="time" ]
		
		#Expand this into an n-dimensional table where each species is 
		#considered jointly with every other (for interactions?)
		#I think this is the right way to expand this. 
		#	1. When you look at s2_temp[1,,], it should be like a matrix of pairwise interactions
		#		per environment 
		#	2. If these were e.g. invasions, then s2_temp[,1,] gives species 1 invading into 
		#		every other species, while s2_temp[,,1] shows 1 as resident and species invading
		#		into 1. 
		s2_temp = array(as.matrix(s1_temp[,2:(nspp+1)] ), dim= c(dim(s1_temp)[1], nspp,nspp) ) 

		#Joint probability distribution: 
		#Now there are two interesting dimensions of information: one in the p(Ni|E), and the other
		#in p(Ni|Nj). Also the simultaneous consideration as p(Ni|E, Nj). 
		#Also, without the actual factorial treatment of both intra and interspecific environment
		#treatments there is an aspect of the information that is missing. 
		s1_joint = prop.table(s2_temp)



	}
}


#2. Make an alphabet out of the 
#Attach all of the output to these variables (Shannon E, Mutual Info, Conditional E)
	SD = NULL
	MI = NULL 
	CE = NULL
	
	nseries = dim(series1)[2] #Assuming each entry is its own time series
	ngens = dim(series1)[1] #Time series length
	
	blocks_k = matrix(0.0, ngens, nseries) # marg of X(k)

	#Need some additional pre-formatting to make this work: 
	#Find number of digits in largest integer: 
	width_A = nchar(trunc(max(series1)))
	k=1 # This should always be 1 for the classic MMI 

	######################
	# p( X(k) )
	######################
	for (f in 1:nseries) { 
		blocks_k[,f] = get_k(series1[,f],k,width_A )
	}
	marg_k =  lapply( (apply(blocks_k,2,table)),  prop.table)

	######################
	# p( X(k),P(k) )
	######################
	joint_kp = c( prop.table(table( as.data.frame(series1) )))

	MMI1 =sum(unlist ( lapply(X_probs,shannon_D2) )) - shannon_D2(X_joint)

an alphabet out of the possible environmental states. For now, just round to
#	2 digits so anything at that level of similarity gets lumped together. 
env_approx = round(unique(mcm1$dynamics.df$env),2)
env_alph = get_k(env_approx, k =1, width_A = 2 )
