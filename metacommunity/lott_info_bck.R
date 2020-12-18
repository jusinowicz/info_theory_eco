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

#=============================================================================
#Numerically solve optimal germination strategies
#=============================================================================
#Germination fraction, in sequence. The endpoints 0 and 1 are special cases 
#which can be avoided. 
H1 = seq(0.05,1,0.05) #Germination fraction.
H1_big= matrix(H1,ngens+1,length(H1),byrow=T) 
Hc = c(matrix("H1",nspp,1) )
#All of the combinations of H1 across species for multispecies competition
Hs_big= eval(parse(text=paste("expand.grid(", paste(unlist(Hc),collapse=","), ")" )))

#Compare these to (uniform) random germination:


#For the average growth rate, rho
env_fit$rho1 = array(0, dim = c(ngens+1, nspp, dim(Hs_big)[1] ) ) #Model 1
env_fit$rho2 = array(0, dim = c(ngens+1, nspp, dim(Hs_big)[1] ) ) #Model 2
env_fit$rho3 = array(0, dim = c(ngens+1, nspp,length(H1) ) ) #Single species

env_fit$m1 = matrix (0, dim(Hs_big)[1],nspp)
env_fit$m2 = matrix (0, dim(Hs_big)[1],nspp)
env_fit$m3 =matrix (0, length(H1),nspp)

#The probability distribution of rho:
breaks = 15
env_fit$pr1 = array(0, dim = c(breaks-1, nspp, dim(Hs_big)[1] ) )
env_fit$pr2 = array(0, dim = c(breaks-1, nspp, dim(Hs_big)[1] ) )
env_fit$pr3 = array(0, dim = c(breaks-1, nspp,length(H1) ) )

#The breaks, which correspond to the rhos/lambdas.
env_fit$br1 = array(0, dim = c(breaks-1, nspp, dim(Hs_big)[1] ) )
env_fit$br2 = array(0, dim = c(breaks-1, nspp, dim(Hs_big)[1] ) )
env_fit$br3 = array(0, dim = c(breaks-1, nspp,length(H1) ) )



for (s in 1:nspp) { 
	
	#Model 1
	for(h in 1:dim(Hs_big)[1]) {

		Hs = as.matrix(unlist(Hs_big[h,]))
		tot_comp = ( env_fit$fr *  env_fit$Ni ) %*% (Hs)
		env_fit$rho1[,s,h] = c(env_fit$sr[s]*(1- Hs[s])   + 
			(1-rowSums( env_fit$sr * env_fit$Ni) ) * env_fit$fr[,s]*(Hs[s])/tot_comp )
		
		#Probability distribution of growth rates
		b_use = seq(min(env_fit$rho1[,s,h]),max(env_fit$rho1[,s,h]), length.out=breaks)
		rho_dist = hist(env_fit$rho1[,s,h],breaks=b_use,plot = FALSE)
		env_fit$pr1[,s,h] = rho_dist$counts/sum(rho_dist$counts)
		env_fit$br1[,s,h] = rho_dist$mids

		#Average growth rate:
		env_fit$m1[h,s] = sum(env_fit$pr1[,s,h]*log(env_fit$br1[,s,h] ) )
	}


	#Model 2
	for(h in 1:dim(Hs_big)[1]) {

		Hs = as.matrix(unlist(Hs_big[h,]))
		tot_comp = ( env_fit$fr *  env_fit$Ni2 ) %*% (Hs)
		env_fit$rho2[,s,h] = c(env_fit$sr[s]*(1- Hs[s])   + env_fit$fr[,s]*(Hs[s])/tot_comp )

		#Probability distribution of growth rates
		b_use = seq(min(env_fit$rho2[,s,h]),max(env_fit$rho2[,s,h]), length.out=breaks)
		rho_dist = hist(env_fit$rho2[,s,h],breaks=b_use,plot = FALSE)
		env_fit$pr2[,s,h] = rho_dist$counts/sum(rho_dist$counts)
		env_fit$br2[,s,h] = rho_dist$mids

		#Average growth rate:
		env_fit$m2[h,s] = sum(env_fit$pr2[,s,h]*log(env_fit$br2[,s,h] ) )

	}

	#Model 3: Single-species model 
	fr_big = matrix(env_fit$fr[,s],ngens+1,length(H1))
	env_fit$rho3[,s,] = log(env_fit$fr[,s]%*%(t(H1))+(1-H1_big)*env_fit$sr[s])

	for(i in 1:length(H1)){ 
		#Probability distribution of growth rates
		b_use = seq(min(env_fit$rho3[,s,i]),max(env_fit$rho3[,s,i]), length.out=breaks)
		rho_dist = hist(env_fit$rho3[,s,i],breaks=b_use,plot = FALSE)
		env_fit$pr3[,s,i] = rho_dist$counts/sum(rho_dist$counts)
		env_fit$br3[,s,i] = rho_dist$mids
	}

	#Average growth rate:
	env_fit$m3[,s] = colMeans(env_fit$rho3[,s,])

}

#3D plot for multispecies competition
data1 = data.frame(x = Hs_big[,1], y = Hs_big[,2], z1=env_fit$m1[,1], z2=env_fit$m1[,2] )
plot_ly() %>% add_trace(data = data1,  x=data1$x, y=data1$y, z=data1$z1, type="mesh3d", intensity =data1$z1  ) 
plot_ly() %>% add_trace(data = data1,  x=data1$x, y=data1$y, z=data1$z1,type="contour" ) 

#Consider intrinsic growth rates jointly: 
# mtx = matrix(NA, nrow=length(unique(data1$x)), ncol=length(unique(data1$y)) )
# mtx[cbind(order(data1$x), order(data1$y))] = data1$z1
mtx1 = data1[order((data1$x)), ]
mtx2 = data1[order((data1$y)), ]
plot(mtx1$z1, mtx2$z2)
data2 =  data.frame(x = mtx1$z1, y = mtx2$z2, z3=mtx1$z1*mtx2$z2)
plot_ly() %>% add_trace(data = data2,  x=data2$x, y=data2$y, z = data2$z3, type="contour" ) 

#=============================================================================
#Instead of solving for optimal germination, find random uniform solution 
#=============================================================================
#Germination fraction, in sequence. The endpoints 0 and 1 are special cases 
#which can be avoided.
nsamp = 1000 
H_runif = matrix( runif(nsamp*nspp), nsamp, nspp) #Germination fraction.


H1_big= matrix(H1,ngens+1,length(H1),byrow=T) 
Hc = c(matrix("H1",nspp,1) )
#All of the combinations of H1 across species for multispecies competition
Hs_big= eval(parse(text=paste("expand.grid(", paste(unlist(Hc),collapse=","), ")" )))

#Compare these to (uniform) random germination:


#For the average growth rate, rho
env_fit$rho_runif1 = array(0, dim = c(ngens+1, nspp, nsamp ) ) #Model 1
env_fit$rho_runif2 = array(0, dim = c(ngens+1, nspp, nsamp ) ) #Model 2
env_fit$rho_runif3 = array(0, dim = c(ngens+1, nspp, nsamp ) ) #Single species

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



for (n in 1:nsamp) { 
	for (s in 1:nspp) { 
		
		Hs = as.matrix(unlist(H_runif[n,]))

		#Model 1
			tot_comp1 = ( env_fit$fr *  env_fit$Ni ) %*% (Hs)
			env_fit$rho_runif1[,s,n] = c(env_fit$sr[s]*(1- Hs[s])   + 
				(1-rowSums( env_fit$sr * env_fit$Ni) ) * env_fit$fr[,s]*(Hs[s])/tot_comp1 )
			
			#Probability distribution of growth rates
			b_use = seq(min(env_fit$rho_runif1[,s,n]),max(env_fit$rho_runif1[,s,n]), length.out=breaks)
			rho_dist = hist(env_fit$rho_runif1[,s,n],breaks=b_use,plot = FALSE)
			env_fit$prunif1[,s,n] = rho_dist$counts/sum(rho_dist$counts)
			env_fit$brunif1[,s,n] = rho_dist$mids

			#Average growth rate:
			env_fit$mr1[n,s] = sum(env_fit$prunif1[,s,n]*log(env_fit$brunif1[,s,n] ) )



		#Model 2

			tot_comp2 = ( env_fit$fr *  env_fit$Ni2 ) %*% (Hs)
			env_fit$rho_runif2[,s,n] = c(env_fit$sr[s]*(1- Hs[s])   + env_fit$fr[,s]*(Hs[s])/tot_comp2 )

			#Probability distribution of growth rates
			b_use = seq(min(env_fit$rho_runif2[,s,n]),max(env_fit$rho_runif2[,s,n]), length.out=breaks)
			rho_dist = hist(env_fit$rho_runif2[,s,n],breaks=b_use,plot = FALSE)
			env_fit$prunif2[,s,n] = rho_dist$counts/sum(rho_dist$counts)
			env_fit$brunif2[,s,n] = rho_dist$mids

			#Average growth rate:
			env_fit$mr2[n,s] = sum(env_fit$prunif2[,s,n]*log(env_fit$brunif2[,s,n] ) )

		#Model 3: Single-species model 
			env_fit$rho_runif3[,s,] = (env_fit$fr[,s]*Hs[s]+(1-Hs[s])*env_fit$sr[s])

			#Probability distribution of growth rates
			b_use = seq(min(env_fit$rho_runif3[,s,n]),max(env_fit$rho_runif3[,s,n]), length.out=breaks)
			rho_dist = hist(env_fit$rho_runif3[,s,n],breaks=b_use,plot = FALSE)
			env_fit$prunif3[,s,n] = rho_dist$counts/sum(rho_dist$counts)
			env_fit$brunif3[,s,n] = rho_dist$mids
			
			#Average growth rate:
			env_fit$mr3[,s] = sum(env_fit$prunif3[,s,n]*log(env_fit$brunif3[,s,n] ) )

	}
}

#3D plot for multispecies competition
data1 = data.frame(x = Hs_big[,1], y = Hs_big[,2], z1=env_fit$m1[,1], z2=env_fit$m1[,2] )
plot_ly() %>% add_trace(data = data1,  x=data1$x, y=data1$y, z=data1$z1, type="mesh3d", intensity =data1$z1  ) 
plot_ly() %>% add_trace(data = data1,  x=data1$x, y=data1$y, z=data1$z1,type="contour" ) 

#Consider intrinsic growth rates jointly: 
# mtx = matrix(NA, nrow=length(unique(data1$x)), ncol=length(unique(data1$y)) )
# mtx[cbind(order(data1$x), order(data1$y))] = data1$z1
mtx1 = data1[order((data1$x)), ]
mtx2 = data1[order((data1$y)), ]
plot(mtx1$z1, mtx2$z2)
data2 =  data.frame(x = mtx1$z1, y = mtx2$z2, z3=mtx1$z1*mtx2$z2)
plot_ly() %>% add_trace(data = data2,  x=data2$x, y=data2$y, z = data2$z3, type="contour" ) 


pl_dm1 = colMeans( (fr_big - env_fit$sr[1])/ (env_fit$fr[,1]%*%(t(H1))+(1-env_fit$sr[1])*H1_big ) )


1-env_fit$sr*colMeans(1/(env_fit$fr) )
#=============================================================================
#Information theory
#=============================================================================


