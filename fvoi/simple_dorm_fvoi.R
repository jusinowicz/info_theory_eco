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
source("./info_theory_functions.R")

#=============================================================================
#Declare variables 
#=============================================================================
ngens = 1000 #Time steps
num_states = 10 #Environmental bins or states
nspp = 2

#Survival rates: 
sr = c(0.8,0.8)
sr = c(0.4,0.4)


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
env_sensed = matrix(1, ngens+1,nspp) #Realized environments from sim

sp_fit_i = matrix(1, ngens+1,nspp)
sp_fit_o = matrix(1, ngens+1,nspp)


sp_act = matrix(1, ngens+1,nspp) #Realized environments from sim


#array( matrix(1, ngens+1,nspp), dim = c(ngens+1,num_states,nspp) ) #Realized fitness from sim
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
#1. 
#env_states = make_env_states(num_states)

#1b.
#env_states[env_states<0.02] = 0 ; env_states = env_states/sum(env_states)

#2. Binomial
env_states = rbinom(ngens,num_states, 0.4)
env_states = hist(env_states,0:(num_states))$counts

#3.Poisson
# env_states = rpois(ngens,num_states)
# num_states = max(env_states)
# env_states = hist(env_states,0:(num_states))$counts

#4.Uniform
# env_states = runif(num_states)

env_states = env_states/sum(env_states)

env = sample(x=(0:(num_states-1)), size=ngens, prob = env_states, replace=T)
#env = make_simple_env(env_states,ngens)
env_prob = prop.table(table(factor(env, levels = 0:(num_states-1))))

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
#fm_method = "constant"


#The conditions for fair/subfair odds are different with this model. Ellner and 
#others have shown that the optimal germination fraction is only <1 when 
#      sr*colMeans(1/fs) > 1
fm = matrix(3,nspp,1) # When this is a constant = num_states, fair odds

fs = matrix(0,num_states,nspp)
# for (s in 1:nspp) { fs[,s] = get_species_fit(probs=env_prob, fcor = fs_cor, fm=fm[s], 
# 	method=fm_method )}

mstates=floor(num_states/2)
fs = get_species_fit_pois(mstates, num_states, nspp,fm )

####Germination fraction
#With gs_cor = 1, the germination fraction is optimal, i.e. fractions 
#match the probability of a state occurring. Decreasing this 
#simulates increasingly poor or mismatched information. cor = 0 
#is no information. 
gs_cor = 0.9999
gs_cor_no = 0
gc = matrix(0.5,nspp,1) #For a constant germination fraction -- matches dormancy model
gs_noi = matrix(0,num_states,nspp) #No info
gs_o = matrix(0,num_states,nspp) #Optimal betting (i.e. proportionate)
gf_method = "variable"

for (s in 1:nspp) { 
	gs_noi[,s] = get_species_fraction(probs = env_prob, gcor = gs_cor, gc = gc[s], 
	method="constant" )
}

#####For the optimum single-species constant rate: 
tsize = 1e4
fr_opt = matrix(1, tsize,nspp)
#fr_opt = array( matrix(1, tsize,nspp), dim = c(tsize,num_states,nspp) ) 
env_opt =  matrix(1, tsize,1)
for (t in 1:tsize){
	fit_tmp = get_fit_one(env_states, fs)
	env_opt[t] = fit_tmp$env_act
	fr_opt[t,] = apply(fit_tmp$sp_fit,2,max)
	# fr_opt[t,,] = fit_tmp$sp_fit
}
	
gs_o =  matrix( c(get_single_opt( fr=fr_opt, nspp=nspp, sr = sr )$opts),num_states,nspp,byrow=T) #Optimal 
#gs_io = matrix(c(get_multi_opt(fr=fr_opt, gs=gs nspp=nspp, sr = sr ) ),num_states,nspp,byrow=T )

####Conditional germination fraction i.e. germination with information
#This function creates a table of conditional probabilities based on the
#G(E|C). Needs acc, which is a number between 0 and 1 describing how accurate
#the cue is on average (1 is perfect accuracy, 0 is none)
gec = get_cp(env_states, acc=c(1,1) )

#Make G(C|E) from G(E|C) to simulate population dynamics later:
gj = gec 
gce = gec
for(s in 1:nspp){
	#Joint probability distribution G(E,C) = G(C,E) = G(E|C) * G(C)
	#gj[,,s] = gec[,,s]*matrix(env_states, num_states,num_states,byrow=T ) 
	#For 0 MI, scramble to Cue: 
	gj[,,s] = matrix(runif(num_states^2), num_states,num_states,byrow=T ) 
	#G(C|E) = G(C,E)/G(E)
	gce[,,s] = gj[,,s]/matrix(rowSums(gj[,,s]),num_states,num_states)
	gce[,,s][!is.finite(gce[,,s])] = 0
}

###########THIS GOT ALL MESSED UP SOMEHOW, NEED TO FIX 
#This is key somehow:
incr=0.01
H1 = seq(0.01,.99,incr)
a1=sr[1]*(1-matrix(H1,length(H1),10))+kronecker(H1,t(fs[,1]))
#a1=matrix(sr*(1-H1), length(sr), dim(fs)[1] )+matrix(H1*fs[,1], length(sr),dim(fs)[1],byrow=T)

g_in_e = H1[apply(a1,2,which.max)]
# g_in_e = matrix(g_in_e,env)
# g_in_e2 = g_in_e
#g_in_e = env_states*g_in_e/(sum(env_states*g_in_e)) 
###########THIS GOT ALL MESSED UP SOMEHOW, NEED TO FIX 



#=============================================================================
#Population dynamics
#=============================================================================	
gs = cbind(env_states,env_states)	
for (t in 1:ngens){
	fit_tmp = get_fit_one(env_states, fs)
	env_act[t] = fit_tmp$env_act #Store the env state
	sp_fit_o[t,] = fs[(env_act[t]+1), ]

	#No information
	#rho_noi[t+1, ] = ( sr*(1-gs_noi[c(which_env)])   + sp_fit * gs_noi[c(which_env)]  ) 
	gs_noi = runif(nspp)
	gnoi_fit[t,] = gs_noi 
	rho_noi[t, ] = ( sr*(1-gnoi_fit[t,])   + sp_fit_o[t,]  * gnoi_fit[t,] ) 
	N_noi[t+1,] = N_noi[t, ] * rho_noi[t, ] 
	N_noi[t+1,][N_noi[t+1,]<0] = 0
	
	#Optimal constant betting
	rho_o[t, ] = ( sr*(1-gs_o[(env_act[t]+1),] )   + sp_fit_o[t,] * gs_o[(env_act[t]+1),]  ) 
	No[t+1,] = No[t, ] * rho_o[t, ] 
	No[t+1,][No[t+1,]<0] = 0

	#With information, conditional germination rates
	#The environment that species sensed (i.e. the cue they got), based on G(C|E)
	sp_fit_tmp = matrix((0:(num_states-1)),num_states,nspp)
	for( s in 1:nspp){ 
		env_sensed[t,s] = sample(x=(0:(num_states-1)), size=1, prob =gce[ (env_act[t]+1),,s], replace=T)
		ec = env_sensed[t,s]
		sp_fit_tmp[,s][sp_fit_tmp[,s]!=ec] = -1 #Identify losers
		sp_fit_tmp[,s][sp_fit_tmp[,s]==ec] = fs[,s][sp_fit_tmp[,s]==ec] #Set winning state to its payout
	}
	sp_fit_i[t,] = sp_fit_tmp[sp_fit_tmp>=0]

	#This version is the most similar to the kelly betting example, but does not 
	#make a lot of ecological sense. In particular, it never makes sense to bet
	#everything on a really crappy year under conditions of sub-fair odds. 
	rho_i[t, ] = ( sr*(1-gec[(env_act[t]+1) , , ][ (env_sensed[t,]+1) ] )   + 
				sp_fit_i[t,] * gec[(env_act[t]+1) , , ][ (env_sensed[t,]+1) ]  )

	# #This version uses the germination rate that matches the sensed environment
	# rho_i[t, ] = ( sr*(1-gs[(env_sensed[t,]+1)] )   + 
	#  			sp_fit_i[t,] * gs[(env_sensed[t,]+1)] )

	# #This version uses the germination rate that matches the sensed environment
	# #and assumes that species are selective as to what donditions they germinate
	# rho_i[t, ] = ( sr*(1-gs[(env_sensed[t,]+1)]*g_in_e[(env_sensed[t,]+1)] )   + 
	#  			sp_fit_i[t,] * gs[(env_sensed[t,]+1)]*g_in_e[(env_sensed[t,]+1)] )

	#This version uses the ideal germination rate that matches the sensed environment
	rho_i[t, ] = ( sr*(1-g_in_e[(env_sensed[t,]+1)] )   + 
	 			sp_fit_i[t,] * g_in_e[(env_sensed[t,]+1)] )


	Ni[t+1,] = Ni[t, ] * rho_i[t, ] 
	Ni[t+1,][Ni[t+1,]<0] = 0
}


#Plot the population growth
nn=1:ngens
plot(log(No[,1]),t="l", ylim = c(0,300))
lines(log(N_noi[,1]), col="red")
lines(log(Ni[,1]),col="blue")

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


####The mutual information between the cue and the environment: 
env_freq = prop.table(table(env_act)) #Environment frequency
sE = shannon_D(env_freq) #Shannon entropy

#For species 1: 
env_act[ngens] = 9
env_act[ngens-1] = 8

c_and_e = prop.table(table( data.frame ( e =env_act, c = env_sensed[,1]) ))  #Joint prob between env and cue
sCgivenE = shannon_CE (c_and_e) #Conditional entropy H(C|E)

#Mutual information: 
mI = sE - sCgivenE 

n2 = 1:200
lines(log(exp(nn*mI) ) ,col="green")
mI_sim = mean(log(rho_i)) - mean(log(rho_o))
lines(log(exp(nn*mI_sim) ) ,col="red")

####Save stuff for figures
save(file ="dm_simp.var",Ni, No, N_noi, rho_noi, rho_o, rho_i, gs_o, gj, gce, gec, 
		 sE, sCgivenE, mI, mI_sim,env_act,env_sensed)
