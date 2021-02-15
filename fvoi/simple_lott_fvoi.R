#=============================================================================
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

#=============================================================================
#Stored values, i.e. population dynamics, information metrics
#=============================================================================
Ni = matrix(1, ngens+1,nspp) #Population
env_act = matrix(1, ngens+1,nspp) #Realized environments from sim

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
#=============================================================================

#With gs_cor = 1, the germination fraction is optimal, i.e. fractions 
#match the probability of a state occurring. Decreasing this 
#simulates increasingly poor or mismatched information. cor = 0 
#is no information. 
gs_cor = 0.9999

##################################
###There are two ways to run this, one of which matches the betting example and 
#the other matches the dormancy model example.
#1. With gf_method = variable and fm_method = either variable or  constant: this 
#	matches betting example
#2. With gf_method = constant and fm_method = variable: this matches the dormancy
#	model.

####Germination fraction
gc = matrix(0.5,nspp,1) #For a constant germination fraction -- matches dormancy model
gs = matrix(0,num_states,nspp)
gf_method = "variable"
for (s in 1:nspp) { gs[,s] = get_species_fraction(probs = env_prob, gcor = gs_cor, gc = gc[s], 
	method=gf_method  )}

####Fitness
#With fs_cor = 1, the fitness matches the optimal germination scheme. 
#Rare events could be weighted to be more valuable instaed by making
#this negative.
fs_cor = 0.999 
#With method = "constant", this is just a constant value per state. 
#With method = "variable," matches based on fs_cor
fm_method = "constant"
fm = matrix(num_states-1,nspp,1) # When this is a constant = num_states, fair odds
#fm = matrix(10,nspp,1) # When this is a constant = num_states, fair odds
fs = matrix(0,num_states,nspp)
for (s in 1:nspp) { fs[,s] = get_species_fit(probs=env_prob, fcor = fs_cor, fm=fm[s], 
	method=fm_method )}

#Simulate annual time-steps
for (t in 1:ngens){
	#Simulate the environment:  
	env_current = apply(env_states, 1, function(x) rbinom(1,100,x) )
	ec = max(env_current)
	env_act[t] = which.max(env_current)

	#Identify species' payoff: 
	sp_fit = matrix(env_current,num_states,nspp)
	sp_fit[sp_fit!=ec] = -1 #Identify losers
	sp_fit[sp_fit==ec] = fs[sp_fit==ec] #Set winning state to its payout

	#New total pop: Betting/germinating proportion * total pop * payout/losses
	#Gi[t+1,] = colSums(gs*Ni[t,]*sp_fit)
	Ni[t+1,] = Ni[t,]*(1+colSums(gs*sp_fit))
	Ni[t+1,][Ni[t+1,]<0] = 0

}

#Plot the population growth
nn=1:ngens
plot(log(Ni[,1]),t="l", ylim = c(0,300))
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
