#=============================================================================
# Functions to go with simple_fvoi 
# These work from the more classic setup of betting on on horse races 
# and build from there. In the classic ecological example of desert annuals, 
# the currency is in seeds, which a population uses to bet on environmental
# states. The payoff is in terms of new seeds (offspring).
# 
#	horses = environmental states. When a horse wins, it is equivalent to 
#			 a particular environment manifesting itself.
#
#	bet    = proportion of population invested in an environmental state. For
#			 example, the number of seeds that germinate. 
#
#   payout = per-capita growth rate. For example, the number of new offspring
#			 (i.e. seeds) produced by germinating individuals. 
#
#=============================================================================

#=============================================================================
# Environment: 
#=============================================================================
#=============================================================================
#	Generate the probabilities for a binomial distribution for a number of 
#	environmental states given by num_states. 
#
#	num_states		Number of discrete environmental states
#=============================================================================
make_env_states = function(num_states=10) {

	counts = rbinom(num_states,10, 0.5)
	probs = counts/sum(counts)
	return(as.matrix(probs) )
}

#=============================================================================
#	Generate a time series where each state has a chance of occurring (i.e.
#	winning) based on generating the max value. 
#
#	num_states		Number of discrete environmental states
#=============================================================================
make_simple_env = function(env_states, ngens = 1000) {

	env_states = as.matrix(env_states)
	num_states = dim(env_states)[1]
	es_fact = factor(1:num_states)

	#Make the temporary environment 
	env_tmp = matrix(0,ngens,num_states) 
	#And this is the final environmental sequence
	env = matrix(0,ngens,1) 

	#Fill each time step from a binomial distribution
	for(n in 1:ngens){ env_tmp[n,] = apply(env_states, 1, function(x) rbinom(1,100,x) ) }

	#Find the winning state
	env = apply(env_tmp,1, which.max )
	env = factor(env,level=es_fact)

	return((env) )
}
#=============================================================================
# Species: 
#=============================================================================
#=============================================================================
#
#=============================================================================
get_species_fraction = function(probs, gcor, gc = 0.5, method = "variable"  ) {

	if (method == "variable") {
		#This is standard code for generating correlated random sequences
		corm = matrix(gcor, nrow = 2, ncol = 2)
	    diag(corm) = 1
	    corm = chol(corm)
	    
	    X2 = runif(length(probs))
	    X = cbind(probs,X2)

	    # induce correlation (does not change X1)
	    new_fraction = X %*% corm
	    #Just take column 2 and renormalize to 1
	    new_fraction = new_fraction[,2]/sum(new_fraction[,2])

	    return( as.matrix( new_fraction) )
	} else if(method == "constant"){
		return(new_fraction = matrix(gc, length(probs),1))
	}
}

#=============================================================================
#
#=============================================================================
get_species_fit = function(probs, fcor, fm, method = "variable" ) {

	if (method == "variable") {
		#This is standard code for generating correlated random sequences
		corm = matrix(fcor, nrow = 2, ncol = 2)
	    diag(corm) = 1
	    corm = chol(corm)
	    
	    X2 = runif(length(probs))
	    X = cbind(probs,X2)

	    # induce correlation (does not change X1)
	    new_fraction = X %*% corm
	    return( (new_fraction[,2]*fm) )

	}else if(method == "constant"){
		return(new_fraction = matrix(fm, length(probs),1))
	}
}


#=============================================================================
# Generate a species fitness distribution that is not correlated with the 
# distribution of environmental states. 
# Use the Poisson distribution, centered on state "mstate", and distributed
# between the min/max environmental state values. 
#=============================================================================
get_species_fit_pois = function(mstates, num_states, nspp,fm ) {
	nsamps = 1e4
	fr_states = matrix(0, nsamps, nspp) 
	new_fraction = matrix(0,num_states,nspp)
	for(s in 1:nspp){
		fr_states[,s] = rpois(nsamps,mstates)
		#renormalize over interval 0, max num_states
		b= num_states-1;a = 0
		fr_states[,s] = (b-a)*(fr_states[,s]-min(fr_states[,s]))/
						(max(fr_states[,s])-min(fr_states[,s])) + 
						a
		new_fraction[,s] = hist(fr_states[,s],breaks = 0:(num_states) )$counts
		new_fraction[,s] = fm[s]* new_fraction[,s]/(max(new_fraction[,s]))
	}

	return(new_fraction)

	# env_states = hist(env_states,0:(num_states))$counts
}

#=============================================================================
#
#=============================================================================
get_fit_one = function(env_states, fs ){
	
	num_states = length(env_states)
	nspp =dim(fs)[2]

	#Simulate the environment:  
	env_current = sample(x=(0:(num_states-1)), size=1, prob =env_states, replace=T)
	ec = max(env_current)
	env_act = ec# which.max(env_current)


	#Identify species' payoff: 
	sp_fit = matrix((0:(num_states-1)),num_states,nspp)
	sp_fit[sp_fit!=ec] = -1 #Identify losers
	sp_fit[sp_fit==ec] = fs[sp_fit==ec] #Set winning state to its payout
	sp_fit[sp_fit<0] = 0 

	fs_env = list( env_act=env_act, sp_fit=sp_fit)

	return(fs_env)
}

#=============================================================================
#get_cp returns the conditional betting proportions of a win based on 
#		having information: b(w|i). Ecologically, this is the "bet" (e.g.
#		germination) on an environment based on a cue: g(e|c). 
#		This function decides the spread of error in information with 
#		the variable acc = [0:1]. When acc = 1, information is perfect and 
#		there is no spread. When acc = 0 there is no information and all 
#		outcomes are equally likely (i.e. uniform). For values in between, 
#		error is generated with exponential decay from the target or true
#		value with more spread as acc -> 0.  
#		acc needs one entry per species 
#=============================================================================

get_cp = function(env_states, acc){ 

	tol = 1e-3
	num_states = length(env_states)
	nspp =length(acc)
	cp_out = array( matrix(0,num_states,num_states), dim = c(num_states,num_states,nspp) )
	xx0=matrix((0:num_states))
 
	for( s in 1:nspp){ 

		#Make an exponential dispersal kernel to show how probability of 
		#error in the conditional probability decays with distance from the
		#env state that matches the cue. With perfect match between cue and 
		#environment (acc =1 ) then there is a single value where e = c. 
		#With no informative cue, this gives a uniform distribution by making
		#a_rr really small (giving kd a large variance). 
		kd=xx0
		a_rr = acc[s]
		if(a_rr == 0 ){a_rr = tol}
		kd = a_rr/2*exp(-a_rr*abs(xx0))
		kd=kd/(sum(kd))

		for (n in 1:num_states){
			cp_tmp = matrix(0,num_states,1)
			cp_er = 1- acc[s] #This is the error, the cp of getting the wrong env
			ll = (num_states+1) - (n) #Distance from the left 
			lr = n   #Distance from the right
			#Get the right side of the kernel 
			cp_tmp[(n+1):num_states] = cp_er/2*kd[2:ll ]
			#Get the left side of the kernel 
			cp_tmp[1:(n-1)] = cp_er/2 * kd[2:lr][(lr-1):1]
			#This is the conditional probability of the match
			cp_tmp[n] = acc[s] 
			if(a_rr == tol ){cp_tmp[n] = mean(cp_tmp); cp_tmp = cp_tmp/sum(cp_tmp) }

			cp_out[,n,s] = cp_tmp[1:num_states]
		}
	}

	return(cp_out)
}

#=============================================================================
#Numerically solve optimal germination strategies for the single-species,
#dormancy model: 
#	Ni[t+1] = Ni[t]( (1-g_i)*s_i + g_i * f_i )
#
#	Ni 			The population vectors (matrix) from the main code
#	env			The vector of environmental states, length must match dim[1] 
#				of Ni. 
#	sr 			Species survival rates, size needs to match dim[2] of Ni.  
#   gw 			An optional weight on the germination rate. This is used to 
#				get the multi-strategy optimum with get_multi_opt
#=============================================================================

get_single_opt= function ( fr, nspp, sr, gw =NULL, incr=0.01) {
	
	ngens = dim(fr)[1]
	#fr = as.matrix(fr)
	nspp = dim(fr)[2]
	if(is.null(gw)) {gw = c(matrix(1,nspp,1))}

	env_fit = NULL
	env_fit$sr =sr
	env_fit$fr = fr
	#Germination fraction, in sequence. The endpoints 0 and 1 are special cases 
	#which can be avoided. 
	H1 = seq(0.01,.99,incr) #Germination fraction.
	Hc = c(matrix("H1",nspp,1) )
	
	#Combinations are independent for singl-species model
	Hs_big = cbind(H1,H1)
	#Hs_big= eval(parse(text=paste("expand.grid(", paste(unlist(Hc),collapse=","), ")" )))

	#For the average growth rate, rho
	env_fit$rho1 = array(1, dim = c(ngens+1, nspp, dim(Hs_big)[1] ) ) #Model 1

	#Make the population time series match rho variables: 
	env_fit$Nj1 = array(0.1, dim = c(ngens+1, nspp, dim(Hs_big)[1] ) )

	#Average of log rho
	env_fit$m1 = matrix (0, dim(Hs_big)[1],nspp)

	#The probability distribution of rho:
	breaks = 15
	env_fit$pr1 = array(0, dim = c(breaks-1, nspp, dim(Hs_big)[1] ) )
	
	#The breaks, which correspond to the rhos/lambdas.
	env_fit$br1 = array(0, dim = c(breaks-1, nspp, dim(Hs_big)[1] ) )

	for(h in 1:dim(Hs_big)[1]) {

		Hs = as.matrix(unlist(Hs_big[h,]))

		#=============================================================================
		#Population dynamics
		#=============================================================================		
		for (n in 1:ngens){
			#Model 3:
			env_fit$rho1[n,,h ] = ( ( env_fit$sr*(1- gw*Hs) )  + 
								env_fit$fr[n,] * gw*Hs) #/(env_fit$fr[n,]*Hs * env_fit$Nj3[n,,h]) )  

			env_fit$Nj1[n+1,,h ] = env_fit$Nj1[n,,h ] * env_fit$rho1[n,,h ] 
		}

		#=============================================================================		
		#Get the optimum
		#=============================================================================		
		env_fit$rho1[,,h] =log(env_fit$rho1[,,h]) 
		env_fit$rho1[,,h][!is.finite(env_fit$rho1[,,h] )] = NA

		for (s in 1:nspp) { 

			#Probability distribution of growth rates
			b_use = seq(min(env_fit$rho1[,s,h],na.rm=T),max(env_fit$rho1[,s,h],na.rm=T), length.out=breaks)
			rho_dist = hist(env_fit$rho1[,s,h],breaks=b_use,plot = FALSE)
			env_fit$pr1[,s,h] = rho_dist$counts/sum(rho_dist$counts)
			env_fit$br1[,s,h] = rho_dist$mids

			#Average log growth rate:
			env_fit$m1[h,s] = sum(env_fit$pr1[,s,h]*(env_fit$br1[,s,h] ) )
		}



	}

	#Which is the max value of the log growth rate in each column? 
	opts=NULL
	opts$opts = Hs_big[ apply(env_fit$m1, 2 ,which.max) ] 
	opts$gr = env_fit$m1[ apply(env_fit$m1, 2 ,which.max) ] 
	return(opts)

}

#=============================================================================
#Numerically solve optimal germination strategies for  single-species
#dormancy model when there are multiple germination strategies. This just 
#applies the single : 
#	Ni[t+1] = Ni[t]( (1-g_i)*s_i + g_i * f_i )
#
#	Ni 			The population vectors (matrix) from the main code
#	env			The vector of environmental states, length must match dim[1] 
#				of Ni. 
#	sr 			Species survival rates, size needs to match dim[2] of Ni.  
#=============================================================================

get_multi_opt= function ( fr, gs, sr, incr=0.01) {
	ngens = dim(fr)[1]
	#fr = as.matrix(fr)
	num_states = dim(gs)[1]
	nspp = dim(gs)[2]

	env_fit = NULL
	env_fit$sr =sr
	env_fit$fr = fr

	gs_io = matrix(0, num_states,nspp) #Optimal germination rates.
	grs = matrix(0, num_states,nspp) #Growth rate of optimum.

	for(g in 1:num_states){
		print( paste("State number:", g) )
		gout = get_single_opt( fr, nspp, sr, gw = gs[g,], incr=0.01)
		gs_io[g,] = gout$opts
		grs[]
	}

	#Germination fraction, in sequence. The endpoints 0 and 1 are special cases 
	#which can be avoided. 
	H1 = seq(0.01,.99,incr) #Germination fraction.
	Hc = c(matrix("H1",nspp,1) )
	
	#Combinations are independent for singl-species model
	Hs_big = cbind(H1,H1)
	#Hs_big= eval(parse(text=paste("expand.grid(", paste(unlist(Hc),collapse=","), ")" )))

	#For the average growth rate, rho
	env_fit$rho1 = array(1, dim = c(ngens+1, nspp, dim(Hs_big)[1] ) ) #Model 1

	#Make the population time series match rho variables: 
	env_fit$Nj1 = array(0.1, dim = c(ngens+1, nspp, dim(Hs_big)[1] ) )

	#Average of log rho
	env_fit$m1 = matrix (0, dim(Hs_big)[1],nspp)

	#The probability distribution of rho:
	breaks = 15
	env_fit$pr1 = array(0, dim = c(breaks-1, nspp, dim(Hs_big)[1] ) )
	
	#The breaks, which correspond to the rhos/lambdas.
	env_fit$br1 = array(0, dim = c(breaks-1, nspp, dim(Hs_big)[1] ) )

	for(h in 1:dim(Hs_big)[1]) {

		Hs = as.matrix(unlist(Hs_big[h,]))

		#=============================================================================
		#Population dynamics
		#=============================================================================		
		for (n in 1:ngens){
			#Model 3:
			env_fit$rho1[n,,h ] = ( ( env_fit$sr*(1- Hs) )  + 
								env_fit$fr[n,] * Hs) #/(env_fit$fr[n,]*Hs * env_fit$Nj3[n,,h]) )  

			env_fit$Nj1[n+1,,h ] = env_fit$Nj1[n,,h ] * env_fit$rho1[n,,h ] 
		}

		#=============================================================================		
		#Get the optimum
		#=============================================================================		
		env_fit$rho1[,,h] =log(env_fit$rho1[,,h]) 
		env_fit$rho1[,,h][!is.finite(env_fit$rho1[,,h] )] = NA

		for (s in 1:nspp) { 

			#Probability distribution of growth rates
			b_use = seq(min(env_fit$rho1[,s,h],na.rm=T),max(env_fit$rho1[,s,h],na.rm=T), length.out=breaks)
			rho_dist = hist(env_fit$rho1[,s,h],breaks=b_use,plot = FALSE)
			env_fit$pr1[,s,h] = rho_dist$counts/sum(rho_dist$counts)
			env_fit$br1[,s,h] = rho_dist$mids

			#Average log growth rate:
			env_fit$m1[h,s] = sum(env_fit$pr1[,s,h]*(env_fit$br1[,s,h] ) )
		}



	}

	#Which is the max value of the log growth rate in each column? 
	opts = Hs_big[ apply(env_fit$m1, 2 ,which.max) ] 

	return(opts)

}


#=============================================================================
#IGNORE THIS ONE FOR NOW -- This corresponds to the bad file
#Numerically solve optimal germination strategies for the single-species,
#dormancy model: 
#	Ni[t+1] = Ni[t]( (1-g_i)*s_i + g_i * f_i )
#
#	Ni 			The population vectors (matrix) from the main code
#	env			The vector of environmental states, length must match dim[1] 
#				of Ni. 
#	sr 			Species survival rates, size needs to match dim[2] of Ni.  
#=============================================================================

get_single_opt_bad = function ( env, nspp, sr, incr=0.05) {
	
	ngens = dim(env)[1]
	env = as.matrix(env)

	env_fit = NULL
	env_fit$sr =sr
	env_fit$fr = env
	#Germination fraction, in sequence. The endpoints 0 and 1 are special cases 
	#which can be avoided. 
	H1 = seq(0.01,.99,incr) #Germination fraction.
	Hc = c(matrix("H1",nspp,1) )
	
	#Combinations are independent for singl-species model
	Hs_big = cbind(H1,H1)
	#Hs_big= eval(parse(text=paste("expand.grid(", paste(unlist(Hc),collapse=","), ")" )))

	#For the average growth rate, rho
	env_fit$rho1 = array(1, dim = c(ngens+1, nspp, dim(Hs_big)[1] ) ) #Model 1

	#Make the population time series match rho variables: 
	env_fit$Nj1 = array(0.1, dim = c(ngens+1, nspp, dim(Hs_big)[1] ) )

	#Average of log rho
	env_fit$m1 = matrix (0, dim(Hs_big)[1],nspp)

	#The probability distribution of rho:
	breaks = 15
	env_fit$pr1 = array(0, dim = c(breaks-1, nspp, dim(Hs_big)[1] ) )
	
	#The breaks, which correspond to the rhos/lambdas.
	env_fit$br1 = array(0, dim = c(breaks-1, nspp, dim(Hs_big)[1] ) )

	for(h in 1:dim(Hs_big)[1]) {

		Hs = as.matrix(unlist(Hs_big[h,]))

		#=============================================================================
		#Population dynamics
		#=============================================================================		
		for (n in 1:ngens){
			#Model 3:
			env_fit$rho1[n,,h ] = ( ( env_fit$sr*(1- Hs) )  + 
								env_fit$fr[n,] * Hs) #/(env_fit$fr[n,]*Hs * env_fit$Nj3[n,,h]) )  

			env_fit$Nj1[n+1,,h ] = env_fit$Nj1[n,,h ] * env_fit$rho1[n,,h ] 
		}

		#=============================================================================		
		#Get the optimum
		#=============================================================================		
		env_fit$rho1[,,h] =log(env_fit$rho1[,,h]) 
		env_fit$rho1[,,h][!is.finite(env_fit$rho1[,,h] )] = NA

		for (s in 1:nspp) { 

			#Probability distribution of growth rates
			b_use = seq(min(env_fit$rho1[,s,h],na.rm=T),max(env_fit$rho1[,s,h],na.rm=T), length.out=breaks)
			rho_dist = hist(env_fit$rho1[,s,h],breaks=b_use,plot = FALSE)
			env_fit$pr1[,s,h] = rho_dist$counts/sum(rho_dist$counts)
			env_fit$br1[,s,h] = rho_dist$mids

			#Average log growth rate:
			env_fit$m1[h,s] = sum(env_fit$pr1[,s,h]*(env_fit$br1[,s,h] ) )
		}



	}

	#Which is the max value of the log growth rate in each column? 
	opts = Hs_big[ apply(env_fit$m1, 2 ,which.max) ] 

	return(opts)

}


#=============================================================================
# Functions to go with lott_info and lott_info_inv. 
#
#	Making species fitness
#	Making species germination (cue)
#	Making environment
#=============================================================================

#=============================================================================
#	Ennvironment: 
#	1. runif1 		random, uniform interval
#	2. rnorm1		random, Gaussian variance
#	3. urand_each,	this considers the optimum of each species' environement
#		nrand_each	and attempts to create an environmental time series with 
#					mode for each species. This requires the additonal 
#					"mweights" to specify the relative frequency of each mode
#	   mweights		weighting for rand_each				
#=============================================================================
get_env = function (env_fit, method = "runif1" ){   
	
	nspp = dim(env_fit$Ni)[2]
	ngens = dim(env_fit$Ni)[1]

	if( method == "runif1") { 
		min1 = 0
		max1 = 1
		if(!is.null(env_fit$min_max) ) { 
			min1 = min(env_fit$min_max)
			max1 = max(env_fit$min_max) 
		}
		env_fit$env = runif(ngens, min = min1, max=max1 ) 
	}

	if( method == "rnorm1") {

		m_use = mean(env_fit$opt)
		v_use = max(env_fit$var)
		if(!is.null(env_fit$g_mean)) {m_use = env_fit$g_mean}
		if(!is.null(env_fit$g_var)) {v_use = env_fit$g_var}
		env_fit$env = rnorm(ngens, m_use, v_use )
	}

	if( method == "urand_each") {
		env_tmp = NULL

		for(s in 1:nspp) {
			min1 = 0
			max1 = 1
			weights = (1/nspp)
			if(!is.null(env_fit$min_max) ) { 
				min1 = min(env_fit$min_max[s,])
				max1 = max(env_fit$min_max[s,]) 
			}
			if(!is.null(env_fit$weights) ) { weights = env_fit$weights[s]}

			env_tmp = c(env_tmp, runif(ngens*weights, min=min1, max=max1) )
		}
		
		lc = ngens - length(env_tmp) 
		if(lc > 0 ) {env_tmp = c(env_tmp,matrix(0,mean(env_tmp),1 ) ) }
		if(lc < 0 ) {env_tmp = env_tmp[-(1:lc)] }


	}
	
	if( method == "nrand_each") {
		
		env_tmp = NULL
		
		for(s in 1:nspp) {
			
			weights = (1/nspp)
			m_use = mean(env_fit$opt)
			v_use = max(env_fit$var)
			if(!is.null(env_fit$g_mean)) { m_use = env_fit$g_mean[s]}
			if(!is.null(env_fit$g_var)) {v_use = env_fit$g_var[s]}
			if(!is.null(env_fit$weights) ) {weights = env_fit$weights[s]}

			env_tmp = c(env_tmp, rnorm(ngens*weights, m_use, v_use ) )
		}
		
		lc = ngens - length(env_tmp) 
		if(lc > 0 ) {env_tmp = c(env_tmp,matrix(mean(env_tmp),abs(lc),1) ) }
		if(lc < 0 ) {env_tmp = env_tmp[-(1:lc)] }

	}

return(sample(env_tmp) )

}


#=============================================================================
#	Fitness:  
#	1. no_var		species have a single environmental value
#	2. uni_var		variance around an optimum that is uniform (runif)
#	3. norm_var		Gaussian around the optimum
#=============================================================================


get_fitness = function (env_fit) { 

	nspp = dim(env_fit$Ni)[2]
	ngens = dim(env_fit$Ni)[1]

	fit_tmp = env_fit$Ni

	for (s in 1:nspp) { 

		if( env_fit$method == "nrand_each" | env_fit$method == "rnorm1" ) {

			m_use = mean(env_fit$opt)
			v_use = max(env_fit$var)
			if(!is.null(env_fit$g_mean)) {m_use = env_fit$g_mean[s]}
			if(!is.null(env_fit$g_var)) {v_use = env_fit$g_var[s]}
			fit_tmp[,s] = exp(-0.5* ( (env_fit$env-m_use)/(v_use) )^2 )
		}

		if( env_fit$method == "urand_each" | env_fit$method == "runif" ) {

			min1 = 0
			max1 = 1
			if(!is.null(env_fit$min_max) ) { 
				min1 = min(env_fit$min_max[s,])
				max1 = max(env_fit$min_max[s,]) 
			}


			fit_tmp[,s] = 1* as.numeric( (env_fit$env) >= min1 & (env_fit$env) <= max1 )
		}
		
		#Standardize on the interval 0,1
			fit_tmp[,s] = (fit_tmp[,s]-  min(fit_tmp[,s]))/( max(fit_tmp[,s]) - min(fit_tmp[,s]) )


	}

return(fit_tmp)


}

#=============================================================================
#	Germination:	
#	1. g_corr		define germination cue relative to fitness using a 
#					correlation coefficient. g_corr = 1 is perfect prediction ,
#					0 is no correlation, negative values would be harmful
#	2. g_always		always germinate a fraction of seeds. 
#	3.				in progress
#=============================================================================
get_env_cue = function (env_fit, method = "g_corr" ){   
	nspp = dim(env_fit$Ni)[2]
	ngens = dim(env_fit$Ni)[1]

	if( method == "g_corr") { 

		cue_tmp = env_fit$Ni
		
		for(s in 1:nspp) {

			cuse = env_fit$g_corr[s]
			#This method will make a second vector for each species that 
			#is correlated to the environmental response.
			ct = cbind(env_fit$fr[,s], runif(ngens) )
			cor1 = cor(ct)
			#1. make independent matrixes
			chol1 = solve(chol(cor1))
			ct_new =  ct %*% chol1 #Make new independent matrix 

			#2. apply new correlation
			cor_use = matrix( c(1, cuse,cuse,1),2,2 )
			chol2 = chol(cor_use) #Factor this 
			ct_use = ct_new %*% chol2 

			#3. standardize on the interval 0.001, 0.999
			ct_use[,2] = (0.999 - 0.001)*(ct_use[,2]-  min(ct_use[,2]))/( max(ct_use[,2]) - min(ct_use[,2]) )+0.001

			cue_tmp[,s] = ct_use[,2]
		}
	}

	return(cue_tmp)
}
