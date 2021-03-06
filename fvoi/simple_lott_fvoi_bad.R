#=============================================================================
# IGNORE THIS ONE FOR NOW
# This one ends up being a weird case the way it is implemented. Not sure yet
# how it fits in.
#
# R code to measure the fitness value of information for a simple model of 
# germination and dormancy from Cohen 1967, e.g. Ellner 1997. This builds
# from the basic example of Kelly betting on horses (where horses = environments
# in the ecological example)
#
# Each time-step, some proportion of each species' seeds may germinate, grow
# and reproduce. The ultimate fitness payoff of this decision is based on 
# how well the germinating seeds have predicted their environment. The optimal
# strategy is to germinate a proportion of seeds that matches on to the 
# probability of sucessful reproduction. Unlike in pure Kelly betting, the 
# population retains a proportion of its "wealth" (population) as a seed bank.
#
# The population dynamics are: 
#
# 	Ni[t+1] = Ni[t]( (1-g_i)*s_i + g_i * f_i )
#	
# Where g_i is the germination rate (this is b_i in Cover and Thomas) and
# f_i is the fitness (this is o_i in Cover and Thomas)
#
# The goal in this problem is usually to find the optimal constant value of g_i.
# (In Kelly betting, g_i is a bet spread over multiple possible states each 
# time step). This is typically considered G(E) in info theory studies. 
# However, we suggest that the organism has already used information about 
# its environment if it knows or has evolved optimal g_i. 
#
# Then G(E|C) would correspond to a variable g_i that is perfectly correlated 
# with the environment. 
#
# Then the equation for FVOI is: deltaG(E;C) = G(E|C) - G(E)
# 
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
ngens = 200 #Time steps
num_states = 10 #Environmental bins or states
nspp = 2

#Survival rates: 
sr = c(0.8,0.8)

#=============================================================================
#Stored values, i.e. population dynamics, information metrics
#=============================================================================
Ni = matrix(1, ngens+1,nspp) #Population w/ information
N_noi = matrix(1, ngens+1,nspp) #Population no information
No = matrix(1, ngens+1,nspp) #Population optimal germination

rho_i = matrix(1, ngens+1,nspp) #Growth rate w/ information
rho_noi = matrix(1, ngens+1,nspp) #Growth rate, no information
rho_o = matrix(1, ngens+1,nspp) #Growth rate, optimal germination

env_act = matrix(1, ngens+1,1) #Realized environments from sim
sp_fit = matrix(1, ngens+1,nspp) #Realized fitness from sim
gi_fit = matrix(1, ngens+1,nspp) #Realized fitness from sim
go_fit = matrix(1, ngens+1,nspp) #Realized fitness from sim
gnoi_fit = matrix(1, ngens+1,nspp) #Realized fitness from sim

#=============================================================================
#Make environment
# 	The environment consists of discrete bins or states, each with some 
#	probability p of occuring. These are calculated in 2 steps: 
#	1. Randomly generate a series of ps for the given number of states.
#	2. Generate sequences and determine the winner at each time step based 
#		on the largest value generated. 
#		This creates a distribution that is highly correlated to the one 
#		generated in step 1, but differs slightly. 
#	3. Count the probability of seeing a state from the simulated sequence
#=============================================================================
env_states = make_env_states(num_states)
env = make_simple_env(env_states,ngens)
env_prob = prop.table(table(env))

#=============================================================================
#Make species' responses
# 	Species respond by deciding how much to bet on each state/what proportion
#	to germinate. 
#   Each state is also assigned a value as the payoff/per-capita growth rate
#	(fitness). 
# 	We explore 3 scenarios here:
#	1. No information -- betting proportions drawn from uniform distribution.
#	2. Optimal germination -- solve for the optimal singl-species proportion
#	3. With information -- species respond to a cue that helps them predict
#	   the environmentc (g_i is conditional, g_i(E|C)
#		
#=============================================================================
##################################
###There are two ways to run this, one of which matches the betting example and 
#the other matches the dormancy model example.
#1. With gf_method = variable and fm_method = either variable or  constant: this 
#	matches betting example
#2. With gf_method = constant and fm_method = variable: this matches the dormancy
#	model.

####Fitness
#With fs_cor = 1, the fitness matches the optimal germination scheme. 
#Rare events could be weighted to be more valuable instaed by making
#this negative.
fs_cor = 0.999 
fm_method = "variable"

#The conditions for fair/subfair odds are different with this model. Ellner and 
#others have shown that the optimal germination fraction is only <1 when 
#      sr*colMeans(1/fs) > 1
fm = matrix(5,nspp,1) # When this is a constant = num_states, fair odds

fs = matrix(0,num_states,nspp)
for (s in 1:nspp) { fs[,s] = get_species_fit(probs=env_prob, fcor = fs_cor, fm=fm[s], 
	method=fm_method )}

####Germination fraction
#With gs_cor = 1, the germination fraction is optimal, i.e. fractions 
#match the probability of a state occurring. Decreasing this 
#simulates increasingly poor or mismatched information. cor = 0 
#is no information. 
gs_cor = 0.9999
gs_cor_no = 0
gc = matrix(0.5,nspp,1) #For a constant germination fraction -- matches dormancy model
gs_noi = matrix(0,num_states,nspp) #No info
gs_i = matrix(0,num_states,nspp) #With information
gf_method = "variable"

for (s in 1:nspp) { 
	gs_i[,s] = get_species_fraction(probs = env_prob, gcor = gs_cor,  
	method=gf_method  )

	gs_noi[,s] = get_species_fraction(probs = env_prob, gcor = gs_cor, gc = gc[s], 
	method="constant" )
}

#For the optimum single-species constant rate: 
tsize = 1e4
fr_opt =  matrix(1, tsize,nspp)
env_opt =  matrix(1, tsize,1)
for (t in 1:tsize){
	fit_tmp = get_fit_one(env_states, fs)
	env_opt[t] = fit_tmp$env_act
	fr_opt[t,] = apply( fit_tmp$sp_fit,2,max)}

gs_o =  matrix( c(get_single_opt( env=fr_opt, nspp=nspp, sr = sr )),num_states,nspp,byrow=T) #Optimal 


#=============================================================================
#Population dynamics
#=============================================================================		
for (t in 1:ngens){
	fit_tmp = get_fit_one(env_states, fs)
	env_act[t] = fit_tmp$env_act #Store the env state
	which_env = apply(fit_tmp$sp_fit, 2, which.max)#[1:nspp] #Which environment was it
	sp_fit[t,] = apply(fit_tmp$sp_fit,2,max) #Get the fitness for this env
	gi_fit[t,] = gs_i[c(which_env)] 
	go_fit[t,] = gs_o[c(which_env)] 

	#Population growth rates:
	#Information
	rho_i[t, ] = ( sr*(1-gi_fit[t,] )   + sp_fit[t,] * gi_fit[t,]  ) 
	Ni[t+1,] = Ni[t, ] * rho_i[t, ] 
	Ni[t+1,][Ni[t+1,]<0] = 0

	#No information
	#rho_noi[t+1, ] = ( sr*(1-gs_noi[c(which_env)])   + sp_fit * gs_noi[c(which_env)]  ) 
	gs_noi = runif(nspp)
	gnoi_fit[t,] = gs_noi 
	rho_noi[t, ] = ( sr*(1-gnoi_fit[t,])   + sp_fit[t,] * gnoi_fit[t,] ) 
	N_noi[t+1,] = N_noi[t, ] * rho_noi[t, ] 
	N_noi[t+1,][N_noi[t+1,]<0] = 0

	#Optimal 
	rho_o[t, ] = ( sr*(1-go_fit[t,])   + sp_fit[t,] * go_fit[t,]  ) 
	No[t+1,] = No[t, ] * rho_o[t, ] 
	No[t+1,][No[t+1,]<0] = 0


}


#Plot the population growth
nn=1:ngens
plot(log(Ni[,1]),t="l", ylim = c(0,300))
lines(log(N_noi[,1]), col="red")
lines(log(No[,1]),col="blue")

#Theoretical prediction based on optimal germination/betting strategy (gs)
lines(log(2^(nn*sum(env_prob*log2(gs[,1]*fs[,1])))),col="red")
#Theoretical prediction when optimal germination matches actual probs
Wbp = log2(env_prob*fs[,1])
Wbp[!is.finite(Wbp)] = 0
lines(log(2^(nn*sum(env_prob*Wbp))),col="blue" )

#The theoretical population growth rate:
#log-Rate: 
lGr = colSums(matrix(env_prob,num_states,nspp)*log(gs*fs))
#Rate:
Gr = apply( (gs*fs)^matrix(env_prob,num_states,nspp),2,prod)
