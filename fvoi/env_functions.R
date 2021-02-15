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
