#=============================================================================
# R code to create a vesion of lottery model where species have cue to respond
# to environments. 
# 
# 1. Create population dynamics based on lottery model
# 2. Species germinate or not in response to a cue. 
# 3. Species reproduction is determined by how well it matches with 
#    environment following germination. 
#=============================================================================
#=============================================================================
#Load libraries
#=============================================================================
library(tidyverse)
library(RandomFields)
library(vegan)
source("./env_functions.R")

#=============================================================================
#Declare variables 
#=============================================================================
ngens = 1000 #Time steps
nspp = 2 #Species

#=============================================================================
#Stored values, i.e. population dynamics, information metrics
#=============================================================================
Ni = matrix(0.1, ngens+1,nspp)

#=============================================================================
#Make environment and species
#
#	There are several ways to measure the environment, and the species' 
#	response to the environment. Choose combinations from these options:
#
#	Fitness:  
#	1. no variance			species have a single environmental value
#	2. uniform variance		variance around an optimum that is uniform (runif)
#	3. normal variance		Gaussian around the optimum
#
#	Germination:	
#	1. g_corr		define germination cue relative to fitness using a 
#					correlation coefficient. g_corr = 1 is perfect prediction ,
#					0 is no correlation, negative values would be harmful
#	2. g_always		always germinate a fraction of seeds. 
#	3.				in progress
#
#	Ennvironment: 
#	1. runif1 		random, uniform interval
#	2. rnorm1		random, Gaussian variance
#	3. urand_each,	this considers the optimum of each species' environement
#		nrand_each	and attempts to create an environmental time series with 
#					mode for each species. This requires the additonal 
#					"mweights" to specify the relative frequency of each mode.
#					urand assumes uniform random, nrand assumes normal 
#	   mweights		weighting for rand_each				
#=============================================================================

#It seems easier to make the species first, then taylor the environment. 

#Fitness: Distribution of species' optimal fitness environment
#Be sure the optimum and variance are chosen with the environmental
#distribution type in mind. 

####1. Make env_fit, species-level properties:
env_fit = NULL
env_fit$Ni = Ni
env_fit$Ni2 = Ni
env_fit$Ni3 = Ni
env_fit$opt = runif(nspp)
env_fit$var = matrix( 0.1 ,nspp,1) #A generic variance
env_fit$min_max = NULL
env_fit$g_mean = NULL
env_fit$g_var = NULL

###Choose 1: 
#A. If the environment is uniform, min_max can be used to specify the upper and lower
#intervals. If the method is "urand_each", this should be a matrix with 
#1 row for each species
# uwidth = 0.1 
# env_fit$min_max = matrix( c( env_fit$opt+uwidth,
# 						 env_fit$opt-uwidth ), nspp,2) #NULL
# #Choose 1:
# #Choose the way that the environment will be constructed: 
# env_fit$method = "runif1"
# env_fit$method = "urand_each"

#B. If the environment is Gaussian, the mean and variance can be specified. Otherwise 
#the average of all species optima is used, and the maximum variance is used. 
#If the method is "nrand_each," these should be vectors with 1 row per species. 
env_fit$g_mean = env_fit$opt 
env_fit$g_var = matrix( 0.1 ,nspp,1)

#Choose 1:
#Choose the way that the environment will be constructed: 
#env_fit$method = "rnorm1"
env_fit$method = "nrand_each"

####2. Environment: Time series of environment. 
#Get the environment:
env_fit$env = get_env(env_fit, method = env_fit$method)

####3. Fitness: Get species' intrinsic fitness in response to the environment. 
env_fit$fr = get_fitness(env_fit)

####4. Cue: Distribution of species' optimal germination environment 
env_fit$cue_method = "g_corr"
#env_fit$cue_dist = "uniform"
#Define how correlated each species' cue is with the environment:
env_fit$g_corr = runif(nspp, min = 0.98, max=0.999)
env_fit$gr= get_env_cue(env_fit, method = env_fit$cue_method)

#Survival rates
env_fit$sr = c(matrix(0.9,nspp,1)) #rnorm(nspp, 0.1, 0.1)

#Scale the intrinsic fitness: 
env_fit$lambda_r = c(5,5)
#Adding a small amount to remove the 0s makes analysis way easier.
env_fit$fr = env_fit$fr* env_fit$lambda_r+.01 

#####
####	Species population-level parameters.  
####	These are a function of the environment. 
#####
#Annual germination rate: 

#Annual intrinsic fitness

#=============================================================================
#Population dynamics
#=============================================================================
for (n in 1:ngens){
	#Lottery model with germination 
	env_fit$Ni[n+1, ] = env_fit$Ni[n, ]*( ( env_fit$sr*(1- env_fit$gr[n,]) )  + 
						(1-sum( env_fit$sr * env_fit$Ni[n, ]) ) * 
						(env_fit$fr[n,]* env_fit$gr[n,]/
					(sum( env_fit$fr[n,]*  env_fit$gr[n,] * env_fit$Ni[n, ]) ) ) ) 

	#"Unscaled" lottery model -- without explicit competition for space
	env_fit$Ni2[n+1, ] = env_fit$Ni2[n, ]*( ( env_fit$sr*(1- env_fit$gr[n,]) )  + 
						(env_fit$fr[n,]* env_fit$gr[n,]/
					(sum( env_fit$fr[n,]*  env_fit$gr[n,] * env_fit$Ni2[n, ]) ) ) ) 

	env_fit$Ni3[n+1, ] = env_fit$Ni3[n, ] * ( ( env_fit$sr*(1- env_fit$gr[n,]) )  + 
						env_fit$fr[n,] * env_fit$gr[n,] ) 

}

#Plot the expected growth rate as a function of germination fraction,
#find the optimal constant germination fraction.

H1 = seq(0,1,0.05) #Germination fraction
H1_big= matrix(H1,ngens+1,length(H1),byrow=T) 
Hc = c(matrix("H1",nspp,1) )

#All of the combinations of H1 across species
Hs_big= eval(parse(text=paste("expand.grid(", paste(unlist(Hc),collapse=","), ")" )))

env_fit$rho1 = array(0, dim = c(ngens+1, nspp, dim(Hs_big)[1] ) ) #Model 1
env_fit$rho2 = array(0, dim = c(ngens+1, nspp, dim(Hs_big)[1] ) ) #Model 2
env_fit$rho3 = array(0, dim = c(ngens+1, nspp,length(H1) ) ) #Single species

for (s in 1:nspp) { 
	
	#Model 1
	for(h in 1:dim(Hs_big)[1]) {

		Hs = unlist(Hs_big[h,])
		tot_comp = ( env_fit$fr *  env_fit$Ni2 ) %*% t(Hs)
		env_fit$rho2[,s,h] = c(env_fit$sr[s]*(1- Hs[s])   + env_fit$fr[,s]*(Hs[s])/tot_comp )

	}


	#Model 2
	for(h in 1:dim(Hs_big)[1]) {

		Hs = unlist(Hs_big[h,])
		tot_comp = ( env_fit$fr *  env_fit$Ni2 ) %*% t(Hs)
		env_fit$rho2[,s,h] = c(env_fit$sr[s]*(1- Hs[s])   + env_fit$fr[,s]*(Hs[s])/tot_comp )

	}

	#Single-species model 
	fr_big = matrix(env_fit$fr[,s],ngens+1,length(H1))
	env_fit$rho3 = log(env_fit$fr[,s]%*%(t(H1))+(1-H1_big)*env_fit$sr[s])


}

eval(parse(paste("expand.grid(", paste(unlist(a2),collapse=","), ")" )))

pl_H1 = colMeans( )
pl_dm1 = colMeans( (fr_big - env_fit$sr[1])/ (env_fit$fr[,1]%*%(t(H1))+(1-env_fit$sr[1])*H1_big ) )

( ( env_fit$sr[1]*(1- H1_big) )  + 
						(env_fit$fr[,1]%*% (t(H1)))/
					(sum( (env_fit$fr* env_fit$Ni2)%*% t(H1) ) ) ) ) 

1-env_fit$sr*colMeans(1/(env_fit$fr) )
#=============================================================================
#Information theory
#=============================================================================


