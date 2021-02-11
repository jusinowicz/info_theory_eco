#=============================================================================
# R code to measure the fitness value of information  for a simple model of 
# germination, made to match the simple examples of proportionate betting
# and bet-hedging
# 
# The general equation for FVOI is deltaG(E;C) = G(E|C) - G(E), where 
# G corresponds to the fitness of a species under a particular setting:
#	G(E|C)	Fitness when species sense the cue (have information).
#	G(E)	Fitness when there is no cue (no information). 
#
# Variability is drawn from a binomial distribution. 
#
# G(E) is calculated as either the optimal rate when there is no cue
# Or, when there is no cue, germination is uniform random form 0 to 1 for all 
#
# G(E|C) is calculated using the actual germination rates. The accuracy of 
# cue/response can be tuned. 
#
# The model is: Ni[t+1] = Ni[t](1+ sum( g_i * f_i )
#	Where g_i is the germination rate (this is b_i in Cover and Thomas) and
#	f_i is the fitness (this is o_i in Cover and Thomas)
#
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
Gi = matrix(1, ngens+1,nspp) #Growth rate


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
gs = matrix(0,num_states,nspp)
for (s in 1:nspp) { gs[,s] = get_species_fraction(env_prob, gs_cor )}
#for (s in 1:nspp) { gs[,s] = env_prob}

#With fs_cor = 1, the fitness matches the optimal germination scheme. 
#Rare events could be weighted to be more valuable instaed by making
#this negative.
#With method = "constant", this is just a constant value per state. 
#With method = "variable," matches based on fs_cor
fs_cor = 0.999 
fm = matrix(num_states-1,nspp,1) # When this is a constant = num_states, fair odds
#fm = matrix(10,nspp,1) # When this is a constant = num_states, fair odds
fs = matrix(0,num_states,nspp)
for (s in 1:nspp) { fs[,s] = get_species_fit(env_prob, fs_cor, fm=fm[s], method="constant" )}
#for (s in 1:nspp) { fs[,s] = get_species_fit(probs = env_prob, fcor = fs_cor, fm=fm[s], method="variable" )}


#Simulate annual time-steps
for (t in 1:ngens){
	#Simulate the environment:  
	env_current = apply(env_states, 1, function(x) rbinom(1,100,x) )
	ec = max(env_current)
	
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
