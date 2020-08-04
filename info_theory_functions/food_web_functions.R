#=============================================================================
# R code to create a simple food web and try to model it using information 
# theoretic tools. 
# 1. Create population dynamics for resource, herbivore, and predator. 
#	 A. Resource is based on a consumer-resource model, with added predators 
#		1. Competition between consumers and resources emerges from consumption
#		2. Parameters at each level can be made a function of temperature. 
#	 B. Resources are stochastic due to environmental fluctuations. 
#	 C. Relative non-linearity allows 2 consumers per Resource
# 2. Change the dynamics by perturbing the resource behavior. 
# 3. Use information theory to track the resulting changes. 
# This left off with not knowing how to represent the decreases in biomass in 
# the compartment model, then finding a chunk of literature about the compart-
# ment models that seems to suggest that models need to be solved first, then
# thinking that the way forward might be to delve into more generic solution 
# techniques first -- i.e. the Lagrangian approach, maximum entropy? 
#=============================================================================
#=============================================================================
# load libraries
#=============================================================================
library(deSolve)
#=============================================================================
# Define a useful function 
#=============================================================================
shifter = function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

#=============================================================================
# food_web_dynamics
# Wrapper funcion for ode in deSolve to set up and run the food web dynamics. 
# Default parameters: 
# spp_list = c(1,1,1)	Number of species at each trophic level
# tend = 1000			Time to simulate over
# delta1 = 0.01			Time step size
# winit = c(3,nspp) 	Initial conditions 
# spp_params			This should be a list with parameters defined for the 
#						the model at each of the 3 trophic levels: 
#						If this is NULL, then it will be generate randomly with
#						each parameter generated as: 
# res_R					This is to specifiy a stochastic forcing function for 
#						the resource growth. The forcing function is lognormal. 
#						It should contain two values:
#						[1] variance and [2] the mean of the undelying normal 
# final 				Only return the last transition (final 2 time steps)? 
#						Do this mainly to reduce memory usage. 			
#=============================================================================

food_web_dynamics = function (spp_list = c(1,1,1), spp_prms = NULL, tend = 1000, 
	delta1 = 0.01, winit =NULL, res_R = NULL, final=FALSE){ 

	###Time to simulate over: 
	times  = seq(from = 0, to = tend, by = delta1)
	tl = length(times)

	###Number of species at each trophic level: 
	nRsp=spp_list[1] #Resource species
	nCsp=spp_list[2] #Consumer species
	nPsp=spp_list[3] #Predator species
	nspp = nRsp+nCsp+nPsp #Total number of species

	###Generic initial conditions: 
	if (is.null(winit) ) {winit = c(matrix(1,nspp,1)) }


	###Parameters of the model:
	if ( length(spp_prms) == 0) {
		spp_prms = NULL
		#Resource: Nearly identical resource dynamics: 
		spp_prms$rR = matrix(rnorm(nRsp,30,0.1), nRsp, 1) #intrinsic growth
		spp_prms$Ki = matrix(rnorm(nRsp,50,0.1), nRsp, 1) #carrying capacity

		#Consumers: 
		spp_prms$rC = matrix(rnorm(nCsp,.8,0.1), nCsp, 1) #intrisic growth
		spp_prms$eFc = matrix(1,nCsp,nRsp) # just make the efficiency for everything 1 for now
		spp_prms$muC = matrix(rnorm(nCsp,0.6,0.1), nCsp, 1) #mortality rates
		#Consumption rates: 
		#Generate a hierarchy where each species predominantly feeds on particular resource. 
		dspp = abs((nCsp - nRsp))
		hier1= seq(1/nRsp, (1-1/nRsp), length=nRsp)
		spp_prms$cC = hier1 
		for( n in 1:(nCsp-dspp)) {
			spp_prms$cC = cbind(spp_prms$cC, shifter(hier1,n))
		}
		spp_prms$cC = as.matrix(spp_prms$cC[1:nRsp,1:nCsp ])

		#Predators: 
		spp_prms$rP = matrix(rnorm(nPsp,0.8,0.1), nPsp, 1) #intrisic growth
		spp_prms$eFp = matrix(1,nPsp,nCsp) # just make the efficiency for everything 1 for now
		spp_prms$muP = matrix(rnorm(nPsp,0.6,0.1), nPsp, 1) #mortality rates
		#Consumption rates: 
		#Generate a hierarchy where each species predominantly feeds on particular resource. 
		dspp = ((nPsp - nCsp))
		if(dspp<0){dspp = 0 }
		hier1= seq(1/nCsp, (1-1/nCsp), length = nCsp)
		spp_prms$cP = hier1
		for( n in 1:(nPsp-dspp-1)) {
			spp_prms$cP = cbind(spp_prms$cP, shifter(hier1,n))
		}
		spp_prms$cP = as.matrix(spp_prms$cP[1:nCsp,1:nPsp])


	}

	rR=(spp_prms$rR); Ki =spp_prms$Ki; rC = spp_prms$rC; eFc = spp_prms$eFc
	muC = spp_prms$muC; cC = spp_prms$cC; rP = spp_prms$rP; eFp = spp_prms$eFp
	muP = spp_prms$muP; cP = spp_prms$cP


	if ( length(res_R)>0){
		a = list( matrix(0,nRsp,1))
		b = list( matrix(0,nCsp,1))

		for( i in 1:nRsp){
			#Make the a variable -- see the documentation for forcings for an example
			amp = res_R[1] #1
			xint = res_R[2] #0
			
			#1. Log normal
			mu = xint
			sd = amp
			location = log(mu^2 / sqrt(sd^2 + mu^2))
			shape = sqrt(log(1 + (sd^2 / mu^2)))
			a[[i]] = approxfun( x = times, y = rlnorm(times, location,shape) , method = "linear", rule = 2) 
			a_t =rlnorm(times, location,shape)
			
			# #2.Normal
			# a[[i]] = approxfun( x = times, y = abs(amp*(rnorm(times) )+xint), method = "linear", rule = 2) 
			# a_t = abs(amp*(rnorm(times) )+xint)

			print( paste("Mean of a(t) = ", mean(a_t),sep="")) 
			print( paste("Var of a(t) = ", var(a_t),sep="")) 
			a_m = mean(a_t)
			vara = var(a_t)
		}

		for( i in 1:nCsp){
			#Make the a variable -- see the documentation for forcings for an example
			amp = res_R[3] #1
			xint = res_R[4] #0
			
			#1. Log normal
			mu = xint
			sd = amp
			location = log(mu^2 / sqrt(sd^2 + mu^2))
			shape = sqrt(log(1 + (sd^2 / mu^2)))
			b[[i]] = approxfun( x = times, y = rlnorm(times, location,shape) , method = "linear", rule = 2) 
			b_t =rlnorm(times, location,shape)
			
			#2.Normal
			# b[[i]] = approxfun( x = times, y = abs(amp*(rnorm(times) )+xint), method = "linear", rule = 2) 
			# b_t = abs(amp*(rnorm(times) )+xint)

			print( paste("Mean of b(t) = ", mean(b_t),sep="")) 
			print( paste("Var of b(t) = ", var(b_t),sep="")) 
			b_m = mean(b_t)
			varb = var(b_t)
		}

	}
	#=============================================================================
	# Define population dynamics
	# This has two different options depending on whether underlying dynamics are
	# stochastic or deterministic. 
	#=============================================================================
	### Function specifying the full dynamics of all trophic levels:

	if ( length(res_R)>0){
		print(res_R)
		#Pass all of these parameters as a list
		parms = list(nspp=nspp, nRsp = nRsp, nCsp = nCsp, nPsp =nPsp,
			rR = spp_prms$rR, Ki =spp_prms$Ki,
			rC = spp_prms$rC, eFc = spp_prms$eFc, muC = spp_prms$muC, cC = spp_prms$cC,
			rP = spp_prms$rP, eFp = spp_prms$eFp, muP = spp_prms$muP, cP = spp_prms$cP,
			a = a
		 )


		food_web = function(times,sp,parms){
				nspp = parms$nspp 
			     R = matrix(sp[1:nRsp], nRsp, 1)
			     C = matrix(sp[(nRsp+1):(nRsp+nCsp)], nCsp, 1)
			     P = matrix(sp[(nRsp+nCsp+1):nspp], nPsp, 1)


				###Resource dynamics: Logistic growth, reduced by consumption
				dR = R
				for( i in 1:nRsp){
					#Logistic - LV consumption
					dR[i] = R[i]*( (rR[i]) * (a[[i]](times)  - R[i]/Ki[i]) - (t(cC[i,])%*%C) )

					#dR[i] = a[[i]](times) + R[i]*( (rR[i]) * (1 - R[i]/Ki[i]) - (t(cC[i,])%*%C))
					
					#Logistic - Saturating consumption
					# dR[i] = a[[i]](times) + R[i]*( (rR[i]) * (1 - R[i]/Ki[i]) - 
					# 	(t(cC[i,])%*%C)^2/(rC[i]^2+(t(cC[i,])%*%C)^2 )  
					# 	)

					#Logistic - Saturating consumption
					# dR[i] = a[[i]](times) + R[i]*( (rR[i]) * (1 - R[i]/Ki[i]) - 
					# 	(t(cC[i,])%*%C)/(rC[i]+(t(cC[i,])%*%C) )  
					# 	)

				}

				###Consumer dynamics
				dC = C 
				for( i in 1:nCsp){
					#LV consumption
					dC[i] = C[i] * ( rC[i] *(b[[i]](times)*eFc[i]*cC[,i])%*%R -(t(cP[i,])%*%P)- muC[i] )
					#dC[i] = C[i] * ( rC[i] *(b[[i]](times)*eFc[i]*cC[,i])%*%R -(t(cP[i,])%*%P)- muC[i] )
					
					#Saturating grazing response.
					# dC[i] = C[i] * ( rC[i] * (
					# 	( (eFc[i]*cC[,i])%*%R)^2/( rC[i]^2 + ((eFc[i]*cC[,i])%*%R)^2  ) ) -
					# 	( (t(cP[i,])%*%P)^2/(rP[i]^2+(t(cP[i,])%*%P)^2 ) )- 
					# 	muC[i] )

					# dC[i] = C[i] * ( rC[i] * (
					# 	( (eFc[i]*cC[,i])%*%R)/( rC[i] + ((eFc[i]*cC[,i])%*%R)  ) ) -
					# 	( (t(cP[i,])%*%P)/(rP[i]+(t(cP[i,])%*%P) ) )- 
					# 	muC[i] )

				}

				###Predator dynamics: LV consumption
				dP = P 
				for( i in 1:nPsp){
					#LV prey consumption
					dP[i] = P[i] * ( rP[i] *(eFp[i]*cP[,i])%*%C - muP[i] )

					# #Saturating consumption
					# dP[i] = P[i] * ( rP[i] *
					# 	( (eFp[i]*cP[,i])%*%C )^2/( rP[i]^2 + ((eFp[i]*cP[,i])%*%C)^2 )  - 
					# 	muP[i] )
					
					#Saturating consumption
					# dP[i] = P[i] * ( rP[i] *
					# 	( (eFp[i]*cP[,i])%*%C )/( rP[i] + ((eFp[i]*cP[,i])%*%C) )  - 
					# 	muP[i] )


				}
			# for( i in 1:nRsp){
			# 	a[[i]] = a[[i]](times) 
			# 	#a = (rR[1]+a(times))
			# }

			return(list(c(dR,dC,dP)))

			}  
	} else {

			#Pass all of these parameters as a list
			parms = list(nspp=nspp, nRsp = nRsp, nCsp = nCsp, nPsp =nPsp,
				rR = spp_prms$rR, Ki =spp_prms$Ki,
				rC = spp_prms$rC, eFc = spp_prms$eFc, muC = spp_prms$muC, cC = spp_prms$cC,
				rP = spp_prms$rP, eFp = spp_prms$eFp, muP = spp_prms$muP, cP = spp_prms$cP
			 )


			food_web = function(times,sp,parms){
				nspp = parms$nspp 
			     R = matrix(sp[1:nRsp], nRsp, 1)
			     C = matrix(sp[(nRsp+1):(nRsp+nCsp)], nCsp, 1)
			     P = matrix(sp[(nRsp+nCsp+1):nspp], nPsp, 1)


				###Resource dynamics: Logistic growth, reduced by consumption
				dR = R
				for( i in 1:nRsp){
					#Logistic - LV consumption
					dR[i] = R[i]*( (rR[i]) * (1 - R[i]/Ki[i]) - (t(cC[i,])%*%C))
					
					# #Logistic - Saturating consumption
					# dR[i] = R[i]*( (rR[i]) * (1 - R[i]/Ki[i]) - 
					# 	(t(cC[i,])%*%C)^2/(rC[i]^2+(t(cC[i,])%*%C)^2 )  
					# 	)

					#Logistic - Saturating consumption
					# dR[i] = R[i]*( (parms$rR[i]) * (1 - R[i]/parms$Ki[i]) - 
					# 	(t(parms$cC[i,])%*%C)/(parms$rC[i]+(t(parms$cC[i,])%*%C) )  
					# 	)

				}

				###Consumer dynamics
				dC = C 
				for( i in 1:nCsp){
					#LV consumption
					dC[i] = C[i] * ( rC[i] *(eFc[i]*cC[,i])%*%R -(t(cP[i,])%*%P)- muC[i] )
					
					#Saturating grazing response.
					# dC[i] = C[i] * ( rC[i] * (
					# 	( (eFc[i]*cC[,i])%*%R)^2/( rC[i]^2 + ((eFc[i]*cC[,i])%*%R)^2  ) ) -
					# 	( (t(cP[i,])%*%P)^2/(rP[i]^2+(t(cP[i,])%*%P)^2 ) )- 
					# 	muC[i] )

					# dC[i] = C[i] * ( parms$rC[i] * (
					# 	( (parms$eFc[i]*parms$cC[,i])%*%R)/( parms$rC[i] + ((parms$eFc[i]*parms$cC[,i])%*%R)  ) ) -
					# 	( (t(parms$cP[i,])%*%P)/(parms$rP[i]+(t(parms$cP[i,])%*%P) ) )- 
					# 	parms$muC[i] )

				}

				###Predator dynamics: LV consumption
				dP = P 
				for( i in 1:nPsp){
					#LV prey consumption
					dP[i] = P[i] * ( rP[i] *(eFp[i]*cP[,i])%*%C - muP[i] )

					# #Saturating consumption
					# dP[i] = P[i] * ( rP[i] *
					# 	( (eFp[i]*cP[,i])%*%C )^2/( rP[i]^2 + ((eFp[i]*cP[,i])%*%C)^2 )  - 
					# 	muP[i] )
					
					#Saturating consumption
					# dP[i] = P[i] * ( parms$rP[i] *
					# 	( (parms$eFp[i]*parms$cP[,i])%*%C )/( parms$rP[i] + ((parms$eFp[i]*parms$cP[,i])%*%C) )  - 
					# 	parms$muP[i] )


				}


			return(list(c(dR,dC,dP)))

			}
	}


	#=============================================================================
	# Run the population models.
	#=============================================================================
	out=NULL
	out_temp = ode(y=winit,times=times,func=food_web,parms=parms)
	#Zero out very small values: 
	out_temp[out_temp<1e-3] = 0

	#Return the whole time series or just the last two time steps? 
	if(final == FALSE){
		out$out= out_temp
	} else { 

		out$out = out_temp[(tl-1):tl,]
	}
	out$spp_prms=parms
	


	return(out)

}

#=============================================================================
# test_eq
# Try to guess whether the system has reached an equilibrium. This just tests
# whether the slope (first derivative) differs significantly from zero, where
# a non-zero slope signifies a non-equilibrium state. It returns a matrix 
# where each entry corresponds to a species in the supplied foodweb with a 
# "1" if that species is not in equilibrium, "0" if it is. 
#
# foodweb 			A variable returned by the function food_web_dynamics
# eq_test			Set how far from the end to test the time series
# t_type			There are two test types: "lm_fit" fits a linear model and
#					tests the slope. "deriv" performs a t-test on the first
#					derivative (or its average). 
#=============================================================================

test_eq = function ( foodweb = NULL, eqtest = NULL,t_type = "deriv") {

	tlast = tlast = dim(foodweb$out)[1] #Time series length

	nRsp = foodweb$spp_prms$nRsp #Numbers of species
	nCsp = foodweb$spp_prms$nCsp
	nPsp = foodweb$spp_prms$nPsp
	nspp = nRsp+nCsp+nPsp

	eq_tf = matrix(1,(nspp+1),1) #Store whether eq (0) or not (1)
	eq_tf [1] = 0 
	for (w in 2:(nspp+1) ){ 

		tm1 = foodweb$out[(eqtest:tlast),1] #Time
		pop1 = foodweb$out[(eqtest:tlast),w] #Population time series

		if(t_type == "lm_fit") { 
			sp_lm = lm( pop1~tm1)
			test1 = summary(sp_lm)$coefficients[2,4] 
			test_ce = sp_lm$coefficients[2]
			#This happens if a species is extinct
			if(is.na(test1) ){
				test1 = 1 
			} 

			if( test1 > 0.05 ) { eq_tf[w] = 0 }
		}

		if(t_type == "deriv") { 
			
			test1= t.test(diff(pop1))$p.value 
			if(is.na(test1) ){
				test1 = 1 
			} 
			if( test1 > 0.05 ) { eq_tf[w] = 0 }
		}
	
	}

	return(eq_tf)

}	



#=============================================================================
# Information theoretic assessment of the foodweb.
#=============================================================================
#=============================================================================
# This code adapts the approach of  Rutledge, Basore, and Mulholland 1976
# A caution when interpreting some of these terms and functions: the fij 
# matrix at the heart of Rutledge's approach are not joint probabilities, but 
# conditional probabilities. See model_notes1.txt for some more clarification
# on this point. In application, this means be sure about whether you are 
# passing the joint or conditional probabilities to functions. 
#=============================================================================
#=============================================================================
# shannon_D
# The Shannon entropy. Rutledge calls it the "Thoroughput " diversity.  
# It is the distribution of energetic flow through a food web. 
# Takes a frequency distrubtion as input. 
# 
# freq 			A frequency distribution that is standardized to 1
#=============================================================================

shannon_D = function ( freq = freq) {

	sD =  - sum ( freq*log(freq),na.rm=T )
	return(sD)

}
 

#=============================================================================
# get_ce
# Conditional entropy in the Rutldege web model. This is mostly implemented to
# double check the math and have multiple routes to the answer.
# Note: The code is written to give H(X|Y) but for the Rutledge web we actually
# want H(Y|X). Make sure the t(fij) is being passed to this function! 
#
# freq 			A frequency distribution that is standardized to 1
# fij 			The transition rates from the Rutledge web 
#=============================================================================

get_ce = function (fij=fij) {
	ncol1 = dim(fij)[1]

	fij_x = rowSums(fij) #Marginals of x and y
	fij_y = colSums(fij)
	cp = fij/matrix(fij_x,ncol1,ncol1) #Conditional probability table P(X|Y)
	ce = -sum(fij_x*rowSums((cp)*log(cp),na.rm=T)) #H(X|Y) = sum over x{ p(x)*H(X|Y=x)}

	return(ce)
}


#=============================================================================
# get_mI
# Average mutual information in the Rutldege web model. The MI is interpreted 
# as the uncertainty resolved by knowing food web structure.
# When energy flow is based on pop size, then the "fij" are the per food source
# conversion rate to new population (I think this should be scaled by the freq of  
# the population i.e. e*c*pop_freq vs. e*c only). 
# Note: in this definition, fij is a conditional probability! 
#
# freq 			A frequency distribution that is standardized to 1
# fij 			The transition rates from the Rutledge web 
#=============================================================================

get_mI = function (freq=freq, fij=fij) {
	ncol1 = length(freq)

	mI = NULL
	mI$mean = 0
	mI$per = matrix(0,dim(fij)[1],dim(fij)[2]) #How much does each link contribute? 

	Pi = fij%*%freq #Biomass proportions at next time
	diag(fij) = 0*diag(fij)
	for (k in 1:ncol1){ 
		for(j in 1:ncol1){
			mI$per[k,j] = freq[k]*fij[k,j] *log(fij[k,j]/(Pi[j]) )
			#print(paste ("first ",freq[k]*fij[k,j],"     then    ",fij[k,j]/(Pi[j]) ))
			# if(is.na(mI$per[k,j])){mI$per[k,j]=0 }
			# if(is.infinite(mI$per[k,j])){mI$per[k,j]=0 } 
		}
	}
	
	mI$mean = sum(rowSums(mI$per,na.rm=T),na.rm=T)
	return(mI)
}

#=============================================================================
# get_mI2
# Average mutual information in the Rutldege web model. The MI is interpreted 
# as the uncertainty resolved by knowing food web structure.
# When energy flow is based on pop size, then the "fij" are the per food source
# conversion rate to new population (I think this should be scaled by the freq of  
# the population i.e. e*c*pop_freq vs. e*c only). This is a second implementation
# based on the standard definition of the MI 
#
# fij 			The transition rates from the Rutledge web 
#=============================================================================

get_mI2 = function (fij=fij) {
	ncol1 = dim(fij)[1]

	fij_x = rowSums(fij) #Marginals of x and y
	fij_y = colSums(fij)
	
	mI = NULL

	pxpy=t(t(fij_x))%*%fij_y #This is p(y)*p(x)
	itmp=(fij*log(fij/pxpy)); itmp[is.na(itmp)]=0
	mI$mean = sum(rowSums(itmp, na.rm=T),na.rm=T) #Mutual information
	mI$per = itmp 

	return(mI)
}


#=============================================================================
# rutledge_web
# This function takes the ODEs and converts them to a biomass balance matrix and 
# transition matrix. This version creates a compartment for each "event" where 
# biomass is gained or loss. This includes birth, death, and "inefficiency" in 
# the form of the way that biomass consumed translates to new population biomass. 
# Note: This function contains no defaults because it is meant to be run 
# on the output of food_web_dynamics()
#
# spp_list				Number of species at each trophic level
# pop_ts				Population matrix with time as rows, each species as column
# spp_prms				This should be a list with parameters defined for the 
#						the model at each of the 3 trophic levels: 
#						If this is NULL, then it will be generate randomly with
#						each parameter generated as: 
# model_type			The compartment model can be constructed in one of several
#						ways. 
# if_conditional		Switch to determine which mathematical approach to take
#						See the notes below. FALSE by default. 
#=============================================================================
#=============================================================================
# Sub-routines for rutledge_web:
#=============================================================================
# get_full_model
# Includes external energetic input into resource biomass production
#=============================================================================
#=============================================================================
# get_resource_model
# 
#=============================================================================

rutledge_web = function (spp_list = spp_list, pop_ts=pop_ts, spp_prms=spp_prms,
							model_type = "full", if_conditional = FALSE) {

	tuse = dim(pop_ts)[1]
	nRsp=spp_list[1] #Resource species
	nCsp=spp_list[2] #Consumer species
	nPsp=spp_list[3] #Predator species
	nspp = nRsp+nCsp+nPsp #Total number of species

	ncol1 = 1+nRsp+2*nCsp+2*nPsp+1
	nrow1 = ncol1

	#Read off the parameters from the model
	rR=(spp_prms$rR); Ki =spp_prms$Ki; rC = spp_prms$rC; eFc = spp_prms$eFc
	muC = spp_prms$muC; cC = spp_prms$cC; rP = spp_prms$rP; eFp = spp_prms$eFp
	muP = spp_prms$muP; cP = spp_prms$cP


	rweb = NULL

	###Generate the quantities that describe "energy" (biomass?) flow through
	###food web.
	rweb$fij = array(c(matrix(0,ncol1,nrow1),matrix(0,ncol1,nrow1)), dim = c(ncol1,nrow1,(tuse))) 
	rweb$Qi = matrix(0,ncol1,(tuse))
	rweb$fijQi = array(c(matrix(0,ncol1,nrow1),matrix(0,ncol1,nrow1)), dim = c(ncol1,nrow1,(tuse))) 
	rweb$pb = matrix(0,ncol1,(tuse))
	rweb$qb = matrix(0,ncol1,(tuse))
	rweb$DC = array(c(matrix(0,ncol1,nrow1),matrix(0,ncol1,nrow1)), dim = c(ncol1,nrow1,(tuse))) 
	rweb$ee = matrix(0,ncol1,(tuse))


	#Information theoretic quantities 
	rweb$sD = matrix(0,tuse, 1) #Shannon entropy
	rweb$mI_mean = matrix(0,tuse, 1) #Mutual Information
	rweb$mI_per =array(c(matrix(0,ncol1,nrow1),matrix(0,ncol1,nrow1)), dim = c(ncol1,nrow1,(tuse))) 
	rweb$ce = matrix(0,tuse, 1) #Conditional Entropy
	rweb$ce2 = matrix(0,tuse, 1) #Conditional Entropy, calculated 2nd way as a check
	rweb$mI_mean2 = matrix(0,tuse, 1) #Mutual Information, calculated 2nd way as a check
	
	###For now, loop through each time step. Is there a faster way to do this with matrix math? 
	for(n in 1:(tuse-1)) { 

		#Make a block matrix of biomass transfer for each trophic level 

		R1 = t(matrix(pop_ts[n,1:nRsp],nCsp,nRsp,byrow=T))
		C1r = t(matrix(pop_ts[n,(nRsp+1):(nRsp+nCsp)],nRsp,nCsp,byrow=T))
		C1p = t(matrix(pop_ts[n,(nRsp+1):(nRsp+nCsp)],nPsp,nCsp,byrow=T))
		P1 = t(matrix(pop_ts[n,(nRsp+nCsp+1):(nspp)],nCsp,nPsp,byrow=T))
		frR = t(matrix(rR,nCsp,nRsp,byrow=T))
		frC = t(matrix(rC,nCsp,nRsp))
		frP = t(matrix(rP,nPsp,nCsp))
		fmuC = t(matrix(muC,nCsp,nRsp))
		fmuP = t(matrix(muP,nPsp,nCsp))


		## Calculate the actual biomass flow between each element. 
		## Productions: 
		#Resource production: 
		#diag(rweb$fijQi[1:nRsp,1:nRsp,n]) = (rR)*R1[,1]
		rweb$fijQi[(1+1:nRsp),(1),n] = (rR)*R1[,1]
		rweb$fijQi[(1),(1),n] = sum((rR)*R1[,1] )
		#Consumer production: 
		rweb$fijQi[(1+(1+nRsp):(nRsp+nCsp)),(1+1:nRsp),n] = t((frC*R1*cC)*t(C1r) )
		#Predator production: 
		rweb$fijQi[(1+(1+nRsp+nCsp):nspp),(1+(1+nRsp):(nRsp+nCsp)),n] = t((frP*C1p*cP)*t(P1))

		##Energetic loss: 
		#Resource to consumer: 
		rweb$fijQi[(1+(nspp+1):(nspp+nCsp)),(1+1:nRsp),n ] = t(( R1*cC)*t(C1r) ) - t((frC*R1*cC)*t(C1r) )
		#Consumer to Predator: 
		rweb$fijQi[(1+(nspp+nCsp+1):(nspp+nCsp+nPsp)),(1+(1+nRsp):(nRsp+nCsp)),n ] = t( (C1p*cP)*t(P1) - (frP*C1p*cP)*t(P1) )

		##Mortality loss:
		#Resource: 
		rweb$fijQi[ncol1,(1+1:nRsp),n ] =(rR)*R1[,1]*(R1[,1]/Ki)
		#Consumer:
		rweb$fijQi[ncol1,(1+(1+nRsp):(nRsp+nCsp)),n ] = t((t(C1r)*fmuC)[1,])
		#Predator
		rweb$fijQi[ncol1,(1+(1+nRsp+nCsp):nspp),n ] = (t(P1)*fmuP)[1,]

		#Now make Qi/Pi
		rweb$Qi[,n] = rowSums( rweb$fijQi[,,n]) #*as.numeric(lower.tri(fijQi[,,n])) )
		
		###Now make fij by standardizing:
		#Note, there are 2 different ways to do this! Each way utilizes a slightly
		#different (but related) set of equations to arrive at the infromation 
		#theoretic quantities. 
		if( if_conditional == TRUE){ 
			
			###1) This is more directly in tune with Rutledge's approach. 
			#	This makes fij a conditional probability matrix.
			rweb$fij[,,n] = rweb$fijQi[,,n]/ matrix(rweb$Qi[,n],ncol1,nrow1,byrow=T)
			rweb$fij[,,n][is.na(rweb$fij[,,n])] = 0
			rweb$fij[,,n][is.infinite(rweb$fij[,,n])] = 0
			#Standardize Qi to be proportions: 
			rweb$Qi[,n] = rweb$Qi[,n]/sum(rweb$Qi[,n])
			#Information theoretic quantities 
			rweb$sD[n] = shannon_D(freq = rweb$Qi[,n])
			mI_temp = get_mI (freq =rweb$Qi[,n], fij=rweb$fij[,,n]  )
			rweb$mI_mean[n] = mI_temp$mean
			rweb$mI_per[,,n] = mI_temp$per
			rweb$ce[n] = rweb$sD[n] - rweb$mI_mean[n]

			
		} else {
			
			###2) This is a more general info theory approach.
			#	This makes fij a joint probability matrix. 
			rweb$fij[,,n] = rweb$fijQi[,,n]/ matrix(sum(rweb$Qi[,n],na.rm=T),ncol1,nrow1,byrow=T)
			rweb$fij[,,n][is.na(rweb$fij[,,n])] = 0
			rweb$fij[,,n][is.infinite(rweb$fij[,,n])] = 0
			#Standardize Qi to be proportions: 
			rweb$Qi[,n] = rowSums( rweb$fij[,,n]) 
			#Information theoretic quantities
			rweb$sD[n] = shannon_D(freq = rweb$Qi[,n])
			mI_temp = get_mI2 (fij=rweb$fij[,,n])
			rweb$mI_mean[n] = mI_temp$mean
			rweb$mI_per[,,n] = mI_temp$per
			rweb$ce2[n] = get_ce(fij=t(rweb$fij[,,n]))
			rweb$mI_mean2 = rweb$sD[n] - rweb$ce2[n]
		
		}
		
	}
return (rweb)
}

#=============================================================================
#Noise and notes:
#=============================================================================
#Pi should be the same as Qi[,n+1]. Should also work equally well before or after 
#standardizing Qi.
#Pi[,n] = fij[,,n]%*%Qi[,n] 

#This is just an example taken from a lecture I found online to help make 
#the functions work. 
# #a1 is the joint probability distribution, p(x,y)
# #marginals are then: rowSums(a1) = p(x), colSums(a1) = p(y)
# x1 = rowSums(a1)
# y1 = colSums(a1)

# a1l = log2(a1); a1l[is.na(a1l)] = 0 ; a1l[is.infinite(a1l)] = 0
# je = -sum(rowSums(a1*a1l)) #Joint entropy 

# HX=sum(y1*log2(y1)) #H(X)
# a2=a1/matrix(x1,4,4) #P(X|Y)
# ce = -sum(x1*colSums((a2)*log2(a2),na.rm=T)) #H(X|Y) = sum over x{ p(x)*H(X|Y=x)}

# pxpy=t(t(x1))%*%y1 #This is p(y)*p(x)
# itmp=(a1*log2(a1/pxpy));itmp[is.na(itmp)]=0
# mI = sum(rowSums(itmp)) #Mutual information


#=============================================================================
# ode_to_mb 
# Convert the ODE parameters and output to a mass-balanace equation. This is 
# to make comparisons with the ECOPATH models. 
#
# pop_ts[n,]*rweb$pb[1:nspp,n] * rweb$ee[1:nspp,n]  - 
# rowSums( 
# t(matrix(pop_ts[n,]*rweb$qb[1:nspp,n], nspp,nspp) ) * rweb$DC[1:nspp,1:nspp, n] )
# = 0
# 
# spp_list				Number of species at each trophic level
# pop_ts				Population matrix with time as rows, each species as column
# spp_prms				This should be a list with parameters defined for the 
#						the model at each of the 3 trophic levels: 
#						If this is NULL, then it will be generate randomly with
#						each parameter generated as: 
# model_type			The compartment model can be constructed in one of several
#						ways. 
# if_conditional		Switch to determine which mathematical approach to take
#						See the notes below. FALSE by default. 
#=============================================================================

ode_to_mb = function (spp_list = spp_list, pop_ts=pop_ts, spp_prms=spp_prms,
							model_type = "full", if_conditional = FALSE) {

	tuse = dim(pop_ts)[1]
	nRsp=spp_list[1] #Resource species
	nCsp=spp_list[2] #Consumer species
	nPsp=spp_list[3] #Predator species
	nspp = nRsp+nCsp+nPsp #Total number of species

	ncol1 = 1+nRsp+2*nCsp+2*nPsp+1
	nrow1 = ncol1

	#Read off the parameters from the model
	rR=(spp_prms$rR); Ki =spp_prms$Ki; rC = spp_prms$rC; eFc = spp_prms$eFc
	muC = spp_prms$muC; cC = spp_prms$cC; rP = spp_prms$rP; eFp = spp_prms$eFp
	muP = spp_prms$muP; cP = spp_prms$cP

	rweb = NULL

	###Generate the quantities that describe "energy" (biomass?) flow through
	###food web.
	rweb$pb = matrix(0,ncol1,(tuse))
	rweb$qb = matrix(0,ncol1,(tuse))
	rweb$DC = array(c(matrix(0,ncol1,nrow1),matrix(0,ncol1,nrow1)), dim = c(ncol1,nrow1,(tuse))) 
	rweb$ee = matrix(0,ncol1,(tuse))


	###For now, loop through each time step. Is there a faster way to do this with matrix math? 
	for(n in 1:(tuse-1)) { 

		#Make a block matrix of biomass transfer for each trophic level 

		R1 = t(matrix(pop_ts[n,1:nRsp],nCsp,nRsp,byrow=T))
		C1r = t(matrix(pop_ts[n,(nRsp+1):(nRsp+nCsp)],nRsp,nCsp,byrow=T))
		C1p = t(matrix(pop_ts[n,(nRsp+1):(nRsp+nCsp)],nPsp,nCsp,byrow=T))
		P1 = t(matrix(pop_ts[n,(nRsp+nCsp+1):(nspp)],nCsp,nPsp,byrow=T))
		frR = t(matrix(rR,nCsp,nRsp,byrow=T))
		frC = t(matrix(rC,nCsp,nRsp))
		frP = t(matrix(rP,nPsp,nCsp))
		fmuC = t(matrix(muC,nCsp,nRsp))
		fmuP = t(matrix(muP,nPsp,nCsp))

		## Calculate the biomass flow between each element, following the 
		## ECOPATH mass-balance conceptualization. 

		#Per-capita production: P/B: 
		#Resource production: 
		rweb$pb[(1:nRsp),n] = (rR)/R1[,1]
		#Consumer production: 
		rweb$pb[((1+nRsp):(nRsp+nCsp)),n] = colSums(frC*R1*cC)/C1r[,1]
		#Predator production: 
		rweb$pb[((1+nRsp+nCsp):nspp),n] = colSums(frP*C1p*cP)/P1[,1]
		rweb$pb[,n][!is.finite(rweb$pb[,n])] = 0

		#Per-capita consumption Q/B: 
		#Consumers: 
		rweb$qb[((1+nRsp):(nRsp+nCsp)),n] = colSums(t(C1r)*cC)/C1r[,1] #colSums(R1*cC) #/C1r[,1]
		#Predators: 
		rweb$qb[((1+nRsp+nCsp):nspp),n] = colSums(t(P1)*cP)/P1[,1] #colSums(C1p*cP) #/P1[,1]
		rweb$qb[,n][!is.finite(rweb$qb[,n])] = 0
	
		#Predation rates by species (the dietary matrix) DCij
		#Consumer consumption: 
		rweb$DC [(1:nRsp),((1+nRsp):(nRsp+nCsp)),n] =  (cC*t(C1r)) / 
											matrix(colSums(cC*t(C1r)),nRsp,nCsp,byrow=T ) 
		#Predator consumption: 
		rweb$DC [((1+nRsp):(nRsp+nCsp)),((1+nRsp+nCsp):nspp),n] =  (t(P1)*cP) / 
											matrix(colSums(t(P1)*cP),nCsp,nPsp,byrow=T ) 
		rweb$DC[,,n][!is.finite(rweb$DC[,,n])] = 0

		#Ecotrophic efficiency, the proportion of biomass that is exported 
		#or predated upon (1-EE is "other mortality" )	
		#Resource
		rweb$ee[(1:nRsp),n] = (1-R1[,1]/Ki) #/
		#Consumer
		rweb$ee[((1+nRsp):(nRsp+nCsp)),n] = (1-fmuC[1,]/colSums(frC*R1*cC)  )
		#Predator
		rweb$ee[((1+nRsp+nCsp):nspp),n] = (1-fmuP[1,]/colSums(frP*C1p*cP) )
				

	}

	return (rweb)
}



#=============================================================================
# rutledge_web2
# This function takes real food web data in the form of a mass-balance model
# and calculates the information theoretic quantities in rutledge_web. The
# essential inputs include the biomass balance matrix, the rates of production,
# rates of consumption, and a diet matrix -- i.e. the proportion of each species' 
# diet that another species represents. 
#
# biomass				The biomass matrix
# pb 					The production matrix
# qb 					The consumption matrix
# DC					The dietary matrix  
# ee 					The ecotrophic efficiency (related to mortality rates)
#
# if_conditional		Switch to determine which mathematical approach to take
#						See the notes below. FALSE by default. 
#=============================================================================
#=============================================================================
# Sub-routines for rutledge_web:
#=============================================================================
# get_full_model
# Includes external energetic input into resource biomass production
#=============================================================================
#=============================================================================
# get_resource_model
# 
#=============================================================================

rutledge_web_mb = function (biomass = biomass, pb = pb, qb = qb, DC = DC, ee = ee,  
						 if_conditional = FALSE) {

	nspp = length(biomass) #Number of species or groups
	pp_mat = colSums(DC) ==0 #Primary producers will not have any positive entries
	nRsp = sum(pp_mat)
	nCsp = nspp - nRsp

	#####
	ncol1 = 1+nRsp+2*nCsp+1
	#####
	nrow1 = ncol1

	rweb = NULL

	###The key quantities for calculating biomass flow from the mass-balance
	###model 
	### Sorting by the pp_mat places the resource spp at the left of the matrix.
	rweb$biomass = biomass[names(sort(pp_mat,decreasing=T))]
	rweb$pb = pb[names(sort(pp_mat,decreasing=T))]
	rweb$pb[rweb$pb<0] = 1
	rweb$qb = qb[names(sort(pp_mat,decreasing=T))]
	rweb$qb[rweb$qb<0] = 0
	rweb$DC = DC[names(sort(pp_mat,decreasing=T)),names(sort(pp_mat,decreasing=T))]
	rweb$ee = ee[names(sort(pp_mat,decreasing=T))]

	###Convert these to the correspnding quantities that describe "energy" (biomass?) 
	### flow in the Rutledge/Information theoretic formulation. 
	rweb$fij = matrix(0,ncol1,nrow1)
	rweb$Qi = matrix(0,ncol1,1)
	rweb$fijQi = matrix(0,ncol1,nrow1)

	#The final information theoretic quantities 
	rweb$sD = matrix(0,1, 1) #Shannon entropy
	rweb$mI_mean = matrix(0,1, 1) #Mutual Information
	rweb$mI_per =matrix(0,ncol1,nrow1) 
	rweb$ce = matrix(0,1, 1) #Conditional Entropy
	rweb$ce2 = matrix(0,1, 1) #Conditional Entropy, calculated 2nd way as a check
	rweb$mI_mean2 = matrix(0,1, 1) #Mutual Information, calculated 2nd way as a check
	

	## Calculate the actual biomass flow between each element. 
	## Productions: 
	#Resource production: pb = (rR)/R1[,1], so fijQi = pb*biomass^2
	#diag(rweb$fijQi[1:nRsp,1:nRsp,n]) = pb *biomass^2
	rweb$fijQi[(1+1:nRsp),(1)] = rweb$pb[1:nRsp] * rweb$biomass[1:nRsp] ^2
	rweb$fijQi[(1),(1)] = sum(rweb$pb *rweb$biomass^2)
	
	#Consumer and predator production: pb = colSums(frC*R1*cC)/C1r[,1], so colSums(frC*R1*cC) = pb*biomass
	#In the mass balance formulation there is no real way to back out frC to fully account for 
	#the efficiency of transfer between each species. The best we can do is find an average
	#efficiency across all species (effectively treating frC as constant)

	## Consumption rates matrix, from the dietary matrix, Q/B, and biomass
	cC1 = t(matrix(rweb$qb*rweb$biomass, nspp,nspp) ) * rweb$DC 
	#Consumption
	consumption = colSums( cC1*matrix(rweb$biomass,nspp,nspp) )
	#Production 
	production = rweb$biomass^2*rweb$pb
	#The efficiency/rate of conversion of spp biomass, 
	fr = production/consumption 
	fr[!is.finite(fr)] = 1
	###
	#Some of these numbers are weird. The value of fr should not be >1, but seems to be sometimes? 
	#If it is, add this line to adjust: 
	### 
	fr_true = fr > 1
	if (sum(fr_true)>0){ 
		pb_new = 0.9*consumption / rweb$biomass^2
		rweb$pb[fr_true] = pb_new[fr_true]
		## Recalculate thes quantities
		cC1 = t(matrix(rweb$qb*rweb$biomass, nspp,nspp) ) * rweb$DC 
		consumption = colSums( cC1*matrix(rweb$biomass,nspp,nspp) )
		production = rweb$biomass^2*rweb$pb
		fr = production/consumption 

	}

	fr = t(matrix(fr,nspp,nspp) ) 

	#Combine these to get the food web transfer rates
	rweb$fijQi[(1+(1+nRsp):(nspp)),(1+1:nspp) ] = t(cC1* fr * matrix(rweb$biomass, nspp,nspp) )[(nRsp+1):nspp, ]

	## Energetic loss: This is also based on the assumption that fr is an average efficiency
	## that does not vary by prey species.  
	rweb$fijQi[(1+(nspp+1) ):(2*nspp-nRsp+1),(1+1:nspp)]= t( (cC1*matrix(rweb$biomass,nspp,nspp)) *(1-fr) )[(nRsp+1):nspp, ]
	
	##Mortality loss:
	#Resource: 
	rweb$fijQi[ncol1,(1+1:nspp)] = rweb$biomass^2*rweb$pb*(1- rweb$ee)
	rweb$fijQi[!is.finite(rweb$fijQi)] = 0
	#Now make Qi/Pi

	rweb$Qi = rowSums( rweb$fijQi) #*as.numeric(lower.tri(fijQi[,,n])) )
	
	
	if( if_conditional == TRUE){ 
		
		###1) This is more directly in tune with Rutledge's approach. 
		#	This makes fij a conditional probability matrix.
		rweb$fij= rweb$fijQi/ matrix(rweb$Qi,ncol1,nrow1,byrow=T)
		rweb$fij[is.na(rweb$fij)] = 0
		rweb$fij[is.infinite(rweb$fij)] = 0
		#Standardize Qi to be proportions: 
		rweb$Qi = rweb$Qi/sum(rweb$Qi)
		#Information theoretic quantities 
		rweb$sD = shannon_D(freq = rweb$Qi)
		mI_temp = get_mI (freq =rweb$Qi, fij=rweb$fij  )
		rweb$mI_mean = mI_temp$mean
		rweb$mI_per= mI_temp$per
		rweb$ce = rweb$sD - rweb$mI_mean

		
	} else {
		
		###2) This is a more general info theory approach.
		#	This makes fij a joint probability matrix. 
		rweb$fij = rweb$fijQi/ matrix(sum(rweb$Qi,na.rm=T),ncol1,nrow1,byrow=T)
		rweb$fij[is.na(rweb$fij)] = 0
		rweb$fij[is.infinite(rweb$fij)] = 0
		#Standardize Qi to be proportions: 
		rweb$Qi = rowSums( rweb$fij) 
		#Information theoretic quantities
		rweb$sD = shannon_D(freq = rweb$Qi)
		mI_temp = get_mI2 (fij=rweb$fij)
		rweb$mI_mean= mI_temp$mean
		rweb$mI_per= mI_temp$per
		rweb$ce2 = get_ce(fij=t(rweb$fij))
		rweb$mI_mean2 = rweb$sD - rweb$ce2
	
	}
	
	
return (rweb)
}

