#=============================================================================
# R code to measure the fitness value of information  for a 
# a vesion of the lottery model where species have a cue to respond
# to variable environments that corresponds to the fitness value over the range
# of possible environments.
# 
# The general equation is deltaG(E;C) = G(E|C) - G(E), where 
# G corresponds to the fitness of a species under a particular setting:
#	G(E|C)	Fitness when species sense the cue.
#	G(E)	FItness when there is no cue. 
#
# This code produces fake species whith a distribution of traits for 
# germination and fitness (two separate parameters), an environment
# that varies, and species germination rates and fitness based on mapping 
# traits to the environment in every generation. 
#
# G(E) is calculated as: When there is no cue, germination is uniform random form 0 to 1 for all 
#   species simultaneously. 
# 	Note, this is different from the classic approach which assumes that G(E) 
#	is an optimal bet-hedging strategy for germination that optimizes fitness
#	in the absence of a cue. However, it is our argument that the ESS for 
# 	this strategy is already making use of information -- i.e. that 
#	selection for an ESS is use of information, 
#   since natural selection reflects an information exchange between 
#	the environment and the trait distribution in a population (as Steve 
# 	Frank might say). 
#
# G(E|C) is calculated using the actual germination rates. The accuracy of 
# cue/response can be tuned. 
#
# This  version includes a loop across competition intensity, as controlled by
# the overlap in optimnal habitat. 
#
#=============================================================================
#=============================================================================
#Load libraries
#=============================================================================
library(tidyverse)
library(RandomFields)
library(vegan)
#source("./env_functions.R")
source("./../fvoi/env_functions.R")

#=============================================================================
#Declare variables 
#=============================================================================
ngens = 2000 #Time steps
nspp = 2 #Species

#=============================================================================
#Stored values, i.e. population dynamics, information metrics
#=============================================================================
Ni = matrix(0.1, ngens+1,nspp)

opt_start = 0.15
opt_inc = 0.005
nincs = opt_start/opt_inc
opts = seq(0.0, opt_start, opt_inc)

#Main variables to store information across loops
env_fit = NULL 
#With info
env_fit$mc2_all = matrix (0, nincs,nspp)
env_fit$mc3_all = matrix (0, nincs,nspp)
#Without info
env_fit$mr1_all = matrix (0, nincs,nspp)
env_fit$mr2_all = matrix (0, nincs,nspp)
env_fit$mr3_all =matrix (0, nincs,nspp)

#=============================================================================
#Outer loop
#=============================================================================

for(o in 1:nincs){ 
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

################################################
####1. Make env_fit, species-level properties:
################################################
env_fit$Ni = Ni #Simple population dynamics
env_fit$Ni2 = Ni #Population dynamics of residents only! 
env_fit$Ni3 = Ni #Dormancy model (no competition)
env_fit$opt = c(opt_start - opts[o], opt_start + opts[o] ) #runif(nspp)
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

###Choose 1:
#Choose the way that the environment will be constructed: 
#env_fit$method = "rnorm1"
env_fit$method = "nrand_each"

################################################
####2. Environment: Time series of environment. 
################################################
#Get the environment:
env_fit$env = get_env(env_fit, method = env_fit$method)

################################################
####3. Fitness: Get species' intrinsic fitness in response to the environment. 
################################################
env_fit$fr = get_fitness(env_fit)

################################################
####4. Cue: Distribution of species' optimal germination environment 
################################################
env_fit$cue_method = "g_corr"
#env_fit$cue_dist = "uniform"
#Define how correlated each species' cue is with the environment:
env_fit$g_corr = runif(nspp, min = 0.98, max=0.98)
env_fit$gr= get_env_cue(env_fit, method = env_fit$cue_method)

################################################
####5. Misc
################################################
#Survival rates
env_fit$sr = c(matrix(0.9,nspp,1)) #rnorm(nspp, 0.1, 0.1)

#Scale the intrinsic fitness: 
env_fit$lambda_r = c(10,10)
#Adding a small amount to remove the 0s makes analysis way easier.
env_fit$fr = env_fit$fr* env_fit$lambda_r+.01 

#####################################################
####	Species population-level parameters.  
####	These are a function of the environment. 
################################################

#Annual germination rate: 

#Annual intrinsic fitness

#=============================================================================
#Fitness with cue detection, invasion growth rates 
#=============================================================================

#For the average growth rate, rho
env_fit$rho_c2 = Ni #Model 2
env_fit$rho_c3 = Ni #Single species

#Average of log rho
env_fit$mc1 = matrix (0, 1,nspp)
env_fit$mc2 = matrix (0, 1,nspp)
env_fit$mc3 = matrix (0, 1,nspp)

#The probability distribution of rho:
breaks = 15
env_fit$prc1 = matrix(0, breaks-1, nspp )
env_fit$prc2 = matrix(0, breaks-1, nspp )
env_fit$prc3 =matrix(0, breaks-1, nspp )

#The breaks, which correspond to the rhos/lambdas.
env_fit$brc1 = matrix(0, breaks-1, nspp )
env_fit$brc2 = matrix(0, breaks-1, nspp )
env_fit$brc3 = matrix(0, breaks-1, nspp )

#For plots of invader vs. resident: 
env_fit$Ni[,1] = 0.01; env_fit$Ni[,-1] = 1

for (s in 1:nspp){ 
	for (n in 1:ngens){

		#Model 2: "Unscaled" lottery model for the residents -- without explicit competition for space
		env_fit$Ni2[n+1, -s] = env_fit$Ni2[n,-s ]*( env_fit$sr[-s]*(1- env_fit$gr[n,-s])  + 
							 env_fit$fr[n,-s]* env_fit$gr[n,-s]/
						(1+sum( env_fit$fr[n,-s]*  env_fit$gr[n,-s] * env_fit$Ni2[n,-s ]) ) )

		#IGR
		env_fit$rho_c2[n,s] = ( ( env_fit$sr[s]*(1- env_fit$gr[n,s]) )  + 
							(env_fit$fr[n,s]* env_fit$gr[n,s]/
						(1+sum( env_fit$fr[n,-s]*  env_fit$gr[n,-s] * env_fit$Ni2[n,-s ]) ) ) ) 

		if (s == 1){ 
			#Model 1: "Unscaled" lottery model for all species
			env_fit$Ni[n+1, ] =  env_fit$Ni[n, ]*( env_fit$sr*(1- env_fit$gr[n, ])  + 
								 env_fit$fr[n,]* env_fit$gr[n, ]/
							(1+sum( env_fit$fr[n, ]*  env_fit$gr[n, ] * env_fit$Ni[n, ]) ) )
		}


		#Model 3: Single species
		env_fit$rho_c3[n,s ] = ( ( env_fit$sr[s]*(1- env_fit$gr[n,s]) )  + 
							env_fit$fr[n,s] * env_fit$gr[n,s] ) 


		env_fit$Ni3[n+1,s ] = env_fit$Ni3[n, s] * env_fit$rho_c3[n,s ]

	}
}

env_fit$rho_c2 = log(env_fit$rho_c2 ) 
env_fit$rho_c3 = log(env_fit$rho_c3)

env_fit$rho_c2[!is.finite(env_fit$rho_C2 )] = NA
env_fit$rho_c3[!is.finite(env_fit$rho_c3 )] = NA

for (s in 1:nspp) { 


	#Probability distribution of growth rates
	b_use = seq(min(env_fit$rho_c2[,s],na.rm=T),max(env_fit$rho_c2[,s],na.rm=T), length.out=breaks)
	rho_dist = hist(env_fit$rho_c2[,s],breaks=b_use,plot = FALSE)
	env_fit$prc2[,s] = rho_dist$counts/sum(rho_dist$counts)
	env_fit$brc2[,s] = rho_dist$mids

	#Average log  growth rate:
	env_fit$mc2_all[o,s] = sum(env_fit$prc2[,s]*(env_fit$brc2[,s] ) )

	#Probability distribution of growth rates
	b_use = seq(min(env_fit$rho_c3[,s],na.rm=T),max(env_fit$rho_c3[,s],na.rm=T), length.out=breaks)
	rho_dist = hist(env_fit$rho_c3[,s],breaks=b_use,plot = FALSE)
	env_fit$prc3[,s] = rho_dist$counts/sum(rho_dist$counts)
	env_fit$brc3[,s] = rho_dist$mids

	#Average log growth rate:
	env_fit$mc3_all[o,s] = sum(env_fit$prc3[,s]*(env_fit$brc3[,s] ) )

}


#=============================================================================
#Find random uniform solution -- i.e. germination varies uniformly between 
#0 and 1 every generation. 
#=============================================================================
#Germination fraction, in sequence. The endpoints 0 and 1 are special cases 
#which can be avoided.
nsamp = 100 

#For the average growth rate, rho
env_fit$rho_runif1 = array(1, dim = c(ngens+1, nspp, nsamp ) ) #Model 1
env_fit$rho_runif2 = array(1, dim = c(ngens+1, nspp, nsamp ) ) #Model 2
env_fit$rho_runif3 = array(1, dim = c(ngens+1, nspp, nsamp ) ) #Single species

#Make the population time series match rho variables: 
env_fit$Nj_runif1 = array(0.1, dim = c(ngens+1, nspp, nsamp ) )
env_fit$Nj_runif2 = array(0.1, dim = c(ngens+1, nspp, nsamp ) )
env_fit$Nj_runif3 = array(0.1, dim = c(ngens+1, nspp, nsamp ) )

#Average log growth rate
env_fit$mr1 = matrix (0, nsamp,nspp)
env_fit$mr2 = matrix (0, nsamp,nspp)
env_fit$mr3 =matrix (0, nsamp,nspp)

#The probability distribution of rho:
breaks = 15
env_fit$prunif1 = array(0, dim = c(breaks-1, nspp, nsamp ) )
env_fit$prunif2 = array(0, dim = c(breaks-1, nspp, nsamp ) )
env_fit$prunif3 = array(0, dim = c(breaks-1, nspp, nsamp) )

#The breaks, which correspond to the rhos/lambdas.
env_fit$brunif1 = array(0, dim = c(breaks-1, nspp, nsamp ) )
env_fit$brunif2 = array(0, dim = c(breaks-1, nspp, nsamp ) )
env_fit$brunif3 = array(0, dim = c(breaks-1, nspp, nsamp ) )

env_fit$Nj_runif1[,1,] = 0.01; env_fit$Nj_runif1[,-1,] = 1


for (h in 1:nsamp) { 
	H_runif = matrix( runif(ngens*nspp), ngens, nspp) #Germination fraction.

	#=============================================================================
	#Population dynamics
	#=============================================================================		
	for ( s in 1:nspp) { 
		for (n in 1:ngens){
			Hs = c(as.matrix(unlist(H_runif[n,])))
		

			#Model 2: "Unscaled" lottery model for the residents -- without explicit competition for space
			#Invader species: 
			env_fit$Nj_runif2[n+1,-s,h ] = env_fit$Nj_runif2[n,-s,h ]* ( ( env_fit$sr[-s]*(1- Hs[-s]) )  + 
								(env_fit$fr[n,-s]* Hs[-s]/
							(1+sum( env_fit$fr[n,-s]*  Hs[-s] * env_fit$Nj_runif2[n,-s,h ]) ) ) )

			env_fit$rho_runif2[n,s,h ] = ( ( env_fit$sr[s]*(1- Hs[s]) )  + 
								(env_fit$fr[n,s]* Hs[s]/
							(1+sum( env_fit$fr[n,-s]*  Hs[-s] * env_fit$Nj_runif2[n, -s ,h ]) ) ) )

			if (s == 1){ 
			#Model 1: "Unscaled" lottery model for all species
			env_fit$Nj_runif1[n+1, ,h] = env_fit$Nj_runif1[n,,h]*( env_fit$sr*(1- Hs)  + 
							 env_fit$fr[n, ]* Hs/
						(1+sum( env_fit$fr[n, ]* Hs * env_fit$Nj_runif1[n,,h ]) ) )
			}


			#Model 3: Single species
			env_fit$rho_runif3[n,s,h ] = ( ( env_fit$sr[s]*(1- Hs[s]) )  + 
								env_fit$fr[n,s] * Hs[s]) #/(env_fit$fr[n,]*Hs * env_fit$Nj3[n,,h]) )  

			env_fit$Nj_runif3[n+1,s,h ] = env_fit$Nj_runif3[n,s,h ] * env_fit$rho_runif3[n,s,h ] 


		}
	}

	#=============================================================================
	#Numerically solve optimal germination strategies
	#=============================================================================
	#Note: All of the probablity and histogram approaches below are more accurate 
	#with the log being taken here:
	env_fit$rho_runif1[,,h] =log(env_fit$rho_runif1[,,h]) 
	env_fit$rho_runif2[,,h] = log(env_fit$rho_runif2[,,h]) 
	env_fit$rho_runif3[,,h] = log(env_fit$rho_runif3[,,h])

	env_fit$rho_runif1[,,h][!is.finite(env_fit$rho_runif1[,,h] )] = NA
	env_fit$rho_runif2[,,h][!is.finite(env_fit$rho_runif2[,,h] )] = NA
	env_fit$rho_runif3[,,h][!is.finite(env_fit$rho_runif3[,,h] )] = NA

	for (s in 1:nspp) { 

		# #Probability distribution of growth rates
		# b_use = seq(min(env_fit$rho_runif1[,s,h],na.rm=T),max(env_fit$rho_runif1[,s,h],na.rm=T), length.out=breaks)
		# rho_dist = hist(env_fit$rho_runif1[,s,h],breaks=b_use,plot = FALSE)
		# env_fit$prunif1[,s,h] = rho_dist$counts/sum(rho_dist$counts)
		# env_fit$brunif1[,s,h] = rho_dist$mids

		# #Average log growth rate:
		# env_fit$mr1[h,s] = sum(env_fit$prunif1[,s,h]*(env_fit$brunif1[,s,h] ) )

		#Probability distribution of growth rates
		b_use = seq(min(env_fit$rho_runif2[,s,h],na.rm=T),max(env_fit$rho_runif2[,s,h],na.rm=T), length.out=breaks)
		rho_dist = hist(env_fit$rho_runif2[,s,h],breaks=b_use,plot = FALSE)
		env_fit$prunif2[,s,h] = rho_dist$counts/sum(rho_dist$counts)
		env_fit$brunif2[,s,h] = rho_dist$mids

		#Average log  growth rate:
		env_fit$mr2[h,s] = sum(env_fit$prunif2[,s,h]*(env_fit$brunif2[,s,h] ) )

		#Probability distribution of growth rates
		b_use = seq(min(env_fit$rho_runif3[,s,h],na.rm=T),max(env_fit$rho_runif3[,s,h],na.rm=T), length.out=breaks)
		rho_dist = hist(env_fit$rho_runif3[,s,h],breaks=b_use,plot = FALSE)
		env_fit$prunif3[,s,h] = rho_dist$counts/sum(rho_dist$counts)
		env_fit$brunif3[,s,h] = rho_dist$mids

		env_fit$mr3[h,s] = sum(env_fit$prunif3[,s,h]*(env_fit$brunif3[,s,h] ) )

	}

}

	env_fit$mr2[env_fit$mr2<0] = NA
	env_fit$mr3[env_fit$mr3<0] = NA
	#Average log growth rate:
	env_fit$mr2_all[o,] = colMeans(env_fit$mr2,na.rm=T) 
	env_fit$mr3_all[o,] = colMeans(env_fit$mr3,na.rm=T) 

}
#
#=============================================================================
#The fitness value of information, single species and with competition.
#With random uniform germination. 
#=============================================================================
deltaG1_comp = env_fit$mc2-colMeans(env_fit$mr2,na.rm=T) 
deltaG1_sing = env_fit$mc3-colMeans(env_fit$mr3,na.rm=T) 

#Conditions for bet-hedging: 1/colMeans(env_fit$fr)*env_fit$sr > 1
colMeans(1/env_fit$fr)*env_fit$sr
plot(env_fit$mc2_all[,1]-env_fit$mr2_all[,1],col="red")                                                                  
points(env_fit$mc3_all[,1]-env_fit$mr3_all[,1])   
plot(env_fit$mc2_all[,1],col="red", ylim=c(0,2) )                                                                  
points(env_fit$mr2_all[,1])       
#=============================================================================
# Plot
#=============================================================================
#Change to long format: 
#Get the runif data in the right format: 
llt_r = apply(env_fit$Nj_runif1,c(1,2),mean)
llt_r = env_fit$Nj_runif1[,,6]
llt = data.frame( cbind(1:(ngens+1), env_fit$Ni, llt_r)); names(llt) = c("time", 1:(2*nspp) )
lott_long =llt %>% gather( species, N, 2:(2*nspp+1))

#Each species' time trajectory
ll_sub = subset(lott_long, time < 50)
p0 = ggplot()+ geom_line( data = ll_sub, aes ( x = time, y = N, color = species)  )+ 
ylab("Population")+  scale_y_log10()+
theme(axis.text.x=element_blank(), axis.title.x=element_blank()) #, legend.position = "none") 
p0



#=============================================================================
#Some plots for optimal bet-hedging, Section 2: 
#=============================================================================
library(plotly)

#A basic single-species plot: 
#data1 = data.frame(x = Hs_big[1:360,1], y = Hs_big[1:360,2], z1=diff(env_fit$m3[,1]), z2=diff(env_fit$m3[,2] ) )
data1 = data.frame(x = Hs_big[,1], y = Hs_big[,2], z1=(env_fit$m3[,1]), z2=(env_fit$m3[,2] ) )

#For multi-species competition, 2 species at a time: 
data2 = data.frame(x = Hs_big[,1], y = Hs_big[,2], z1=env_fit$m2[,1], z2=env_fit$m2[,2], z3 = apply(env_fit$m2,1, prod)  )
#data2 = data.frame(x = Hs_big[,1], y = Hs_big[,2], z1=env_fit$m2[,1], z2=env_fit$m2[,2], z3 = apply(env_fit$m2,1, sum)  )

mspp1 =  which(data2$z3 == max(data2$z3,na.rm=T))
mspp2 =  which(round(data2$z3,2) == max(round(data2$z3,2),na.rm=T))

par(mfrow=c(2,1))
plot(data1$x, data1$z1, ylim= c(-0.1,1.5) )
points(data2$x, data2$z1,col="red")
points(data2$x[mspp1], data2$z1[mspp1],col="blue")
points(data2$x, data2$z3, col="green")

plot(data1$y, data1$z2, ylim= c(-0.1,1.5) )
points(data2$y, data2$z2,col="red")
points(data2$y[mspp1], data2$z2[mspp1],col="blue")
points(data2$y, data2$z3, col="green")




#3D plot for single species:
plot_ly() %>% add_trace(data = data1,  x=data1$x, y=data1$y, z=data1$z1, type="mesh3d", intensity =data1$z1  ) 

#
plot_ly() %>% add_trace(data = data1,  x=data1$x, y=data1$y, z=data1$z1,type="contour" ) 


#A 3D plot for multispecies competition
plot_ly() %>% add_trace(data = data2,  x=data2$x, y=data2$y, z=data2$z1, type="mesh3d", intensity =data2$z1  ) 

#
plot_ly() %>% add_trace(data = data2,  x=data2$x, y=data2$y, z=data2$z1,type="contour" ) 

# For the combined metric z3
plot_ly() %>% add_trace(data = data2,  x=data2$x, y=data2$y, z=data2$z3, type="mesh3d", intensity =data2$z3  ) 

#
plot_ly() %>% add_trace(data = data2,  x=data2$x, y=data2$y, z=data2$z3,type="contour" ) 

#Consider intrinsic growth rates jointly: 
# mtx = matrix(NA, nrow=length(unique(data1$x)), ncol=length(unique(data1$y)) )
# mtx[cbind(order(data1$x), order(data1$y))] = data1$z1
# mtx1 = data1[order((data1$x)), ]
# mtx2 = data1[order((data1$y)), ]
# plot(mtx1$z1, mtx2$z2)
# data2 =  data.frame(x = mtx1$z1, y = mtx2$z2, z3=mtx1$z1*mtx2$z2)
# plot_ly() %>% add_trace(data = data2,  x=data2$x, y=data2$y, z = data2$z3, type="contour" ) 

#=============================================================================
#Some plots for random uniform germination, Section 3: 
#=============================================================================
#A basic single-species plot: 
data1b = data.frame(x = H_runif[,1], y = H_runif[,2], z1=env_fit$mr3[,1], z2=env_fit$mr3[,2] )
#For multi-species competition, 2 species at a time: 
data2b = data.frame(x = H_runif[,1], y = H_runif[,2], z1=env_fit$mr2[,1], z2=env_fit$mr2[,2], z3 = apply(env_fit$mr2,1, sum)  )
mspp1b =  which(data2b$z3 == max(data2b$z3,na.rm=T))
mspp2b =  which(round(data2b$z3,2) == max(round(data2b$z3,2),na.rm=T))


par(mfrow=c(2,1))
plot(data1b$x, data1b$z1, ylim= c(-0.1,1.5) )
points(data2b$x, data2b$z1,col="red")
points(data2b$x[mspp1], data2b$z1[mspp1],col="blue")
points(data2b$x, data2b$z3, col="green")

plot(data1b$y, data1b$z2, ylim= c(-0.1,1.5) )
points(data2b$y, data2b$z2,col="red")
points(data2b$y[mspp1], data2b$z2[mspp1],col="blue")
points(data2b$y, data2b$z3, col="green")

#Conditions for bet-hedging: 1/colMeans(env_fit$fr)*env_fit$sr > 1
1/colMeans(env_fit$fr)*env_fit$sr

#3D plot for single species:
plot_ly() %>% add_trace(data = data1b,  x=data1b$x, y=data1b$y, z=data1b$z1, type="mesh3d", intensity =data1b$z1  ) 

#
plot_ly() %>% add_trace(data = data1b,  x=data1b$x, y=data1b$y, z=data1b$z1,type="contour" ) 


#A 3D plot for multispecies competition
plot_ly() %>% add_trace(data = data2b,  x=data2b$x, y=data2b$y, z=data2b$z1, type="mesh3d", intensity =data2b$z1  ) 

#
plot_ly() %>% add_trace(data = data2b,  x=data2b$x, y=data2b$y, z=data2b$z3, type="mesh3d", intensity =data2b$z3  ) 

#
plot_ly() %>% add_trace(data = data1b,  x=data1b$x, y=data1$y, z=data1b$z1,type="contour" ) 

which(data1b$z3 == max(data1b$z3))



#=============================================================================
#Information theory
#=============================================================================


