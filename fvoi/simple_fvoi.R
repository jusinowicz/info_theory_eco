#=============================================================================
# R code to measure the fitness value of information for a simple model of 
# germination, made to match the simple examples of proportionate betting
# and bet-hedging
# 
# This is the simplest example, which doesn't really map on to the biology
# yet, but is a good starting place to test code and conceptual relations. 
#
# Each time-step, all of the population/money is bet on the entire list 
# of possible environmental states. Translated into a population model, 
# it would look something like this: 
#
# 	Ni[t+1] = Ni[t](1+ sum( g_i * f_i )
#	
# Where g_i is the germination rate (this is b_i in Cover and Thomas) and
# f_i is the fitness (this is o_i in Cover and Thomas)
# The sum is over possible environmental states. Thus the product g_i*f_i
# will contain negative values for incorrect guesses/bets. 
# 
#=============================================================================
#=============================================================================
#Load libraries
#=============================================================================
library(tidyverse)
library(RandomFields)
library(vegan)
source("./env_functions.R")
source("../info_theory_functions/info_theory_functions.R")
#=============================================================================
#Declare variables 
#=============================================================================
ngens = 1000 #Time steps
num_states = 10 #Environmental bins or states
nspp = 2

#=============================================================================
#Stored values, i.e. population dynamics, information metrics
#=============================================================================
Ni = matrix(10, ngens+1,nspp) #Population
Ni2 = matrix(10, ngens+1,nspp) #Population
Ni_i = matrix(10, ngens+1,nspp) #Population with information
rhoi = matrix(10, ngens+1,nspp) #Population
rhoi2 = matrix(10, ngens+1,nspp) #Population
rhoi_i = matrix(10, ngens+1,nspp) #Population with information
env_act = matrix(1, ngens+1,1) #Realized environments from sim
env_sensed = matrix(1, ngens+1,nspp) #Realized environments from sim


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
#These are to check numbers and theory: 
# num_states = 3
# env_states = c(0.5,.25,.25)

#env_states = make_env_states(num_states)
env_states = rbinom(ngens,num_states, 0.4)
env_states = hist(env_states,0:(num_states))$counts + 1e-4
env_states = env_states/sum(env_states)

env = sample(x=(0:(num_states-1)), size=ngens, prob =env_states, replace=T)
#env = make_simple_env(env_states,ngens)
env_prob = prop.table(table(factor(env, levels = 0:(num_states-1))))

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
gs_cor = 0.99999

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
fm = matrix(num_states,nspp,1) # When this is a constant = num_states, fair odds
#fm = matrix(10,nspp,1) # When this is a constant = num_states, fair odds
fs = matrix(0,num_states,nspp)
for (s in 1:nspp) { fs[,s] = get_species_fit(probs=env_prob, fcor = fs_cor, fm=fm[s], 
	method=fm_method )}

####Conditional germination fraction i.e. germination with information
#This function creates a table of conditional probabilities based on the
#G(E|C)
gec = get_cp(env_states, acc=c(1,1) )
#Make G(C|E) from G(E|C) to simulate population dynamics later:
gj = gec 
gce = gec
for(s in 1:nspp){
	#Joint probability distribution: G(E|C) * G(C)
	gj[,,s] = gec[,,s]*matrix(gs[,s], num_states,num_states,byrow=T ) 
	#G(C|E) = G(C,E)/G(E)
	gce[,,s] = gj[,,s]/matrix(rowSums(gj[,,s]),num_states,num_states)
}

#Simulate annual time-steps
for (t in 1:ngens){
	#Simulate the environment:  
	#env_current = apply(env_states, 1, function(x) rbinom(1,100,x) )
	env_current = sample(x=(0:(num_states-1)), size=1, prob =env_states, replace=T)
	ec = max(env_current)
	env_act[t] = ec# which.max(env_current)

	#Identify species' payoff: 
	sp_fit = matrix((0:(num_states-1)),num_states,nspp)
	sp_fit[sp_fit!=ec] = -1 #Identify losers
	sp_fit[sp_fit==ec] = fs[sp_fit==ec] #Set winning state to its payout
	sp_fit[sp_fit<0] = 0

	####Without information, optimal
	#New total pop: Betting/germinating proportion * total pop * payout/losses
	Ni[t+1,] = (colSums(matrix(Ni[t,],num_states, nspp,byrow=T)*gs*sp_fit))
	
	#This should be exactly the same: 
	sp_fit2 = sp_fit[sp_fit>0 ]
	gs2 = gs[(ec+1),]
	rhoi2 [t+1, ] = gs2*sp_fit2
	Ni2[t+1,] =  ( Ni2[t,]*gs2*sp_fit2)

	####With information, as determined by conditional probabilities
	#The environment that species sensed (i.e. the cue they got), based on G(C|E)
	sp_fit_i = matrix((0:(num_states-1)),num_states,nspp)
	for( s in 1:nspp){ 
		env_sensed[t,s] = sample(x=(0:(num_states-1)), size=1, prob =gce[ (env_act[t]+1),,s], replace=T)
		ec = env_sensed[t,s]
		sp_fit_i[,s][sp_fit_i[,s]!=ec] = -1 #Identify losers
		sp_fit_i[,s][sp_fit_i[,s]==ec] = fs[,s][sp_fit_i[,s]==ec] #Set winning state to its payout
	}
	sp_fit_i[sp_fit_i<0] = 0 

	rhoi_i[t+1,] = colSums(gec[(env_act[t]+1) , , ][ (env_sensed[t,]+1) ]*
					sp_fit_i)
	Ni_i[t+1,] = (colSums(matrix(Ni_i[t,],num_states, nspp,byrow=T)*
					gec[(env_act[t]+1) , , ][ (env_sensed[t,]+1) ]*
					sp_fit_i))


}

#Plot the population growth
nn=1:ngens
plot(log(Ni[,1]),t="l", ylim = c(0,300))
#Theoretical prediction based on optimal germination/betting strategy (gs)
#lines(log(exp(nn*sum(env_states*log(env_states*fs[,1]) ) ) ),col="red")
#Theoretical prediction when optimal germination matches actual probs
# Wbp = log(env_prob*fs[,1])
# Wbp[!is.finite(Wbp)] = 0
# lines(log(exp(nn*sum(env_prob*Wbp))),col="blue" )

#Add in concditional population growth (growth with cue/information)
lines(log(Ni_i[,1]),t="l", col="blue")

#The theoretical population growth rate:
#log-Rate: 
lGr = colSums(matrix(env_prob,num_states,nspp)*log(gs*fs))
#Rate:
Gr = apply( (gs*fs)^matrix(env_prob,num_states,nspp),2,prod)

####The mutual information between the cue and the environment: 
env_freq = prop.table(table(env_act)) #Environment frequency
sE = shannon_D(env_freq) #Shannon entropy

#For species 1: 
c_and_e = prop.table(table( data.frame ( e =env_act, c = env_sensed[,1]) ))  #Joint prob between env and cue
sCgivenE = shannon_CE (c_and_e) #Conditional entropy H(C|E)

#Mutual information: 
mI = sE - sCgivenE 

n2 = 1:200
lines(log(exp(nn*mI) ) ,col="green")
mI_sim = mean(log(rhoi_i)) - mean(log(rhoi2))
lines(log(exp(nn*mI_sim) ) ,col="red")


