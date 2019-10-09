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
# spp_params			This should be a list with parameters defined for the 
#						the model at each of the 3 trophic levels: 
#						If this is NULL, then it will be generate randomly with
#						each parameter generated as: 
# res_R					This is to specifiy a stochastic forcing function for 
#						the resource growth. The forcing function is lognormal. 
#						It should contain two values:
#						[1] variance and [2] the mean of the undelying normal 
#			
#=============================================================================

food_web_dynamics = function (spp_list = c(1,1,1), spp_prms = NULL, tend = 1000, 
	delta1 = 0.01, res_R = NULL ){ 

	###Time to simulate over: 
	times  = seq(from = 0, to = tend, by = delta1)
	tl = length(times)

	###Number of species at each trophic level: 
	nRsp=spp_list[1] #Resource species
	nCsp=spp_list[2] #Consumer species
	nPsp=spp_list[3] #Predator species
	nspp = nRsp+nCsp+nPsp #Total number of species

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
		
		#Make the a variable -- see the documentation for forcings for an example
		amp = res_R[1] #1
		xint = res_R[2] #0
		a = approxfun( x = times, y = amp*exp(rnorm(times)+xint), method = "linear", rule = 2) 
		a_t = exp(amp*rnorm(times)+xint)
		#a = approxfun( x = times, y = amp*rnorm(times)+xint, method = "linear", rule = 2) 
		#a_t = amp*rnorm(times)+xint
		print( paste("Var in a(t) = ", var(a_t),sep="")) 
		a_m = mean(a_t)
		vara = var(a_t)

	}
	#=============================================================================
	# Define population dynamics
	# This has two different options depending on whether underlying dynamics are
	# stochastic or deterministic. 
	#=============================================================================
	### Function specifying the full dynamics of all trophic levels:

	if ( length(res_R)>0){

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
					dR[i] = R[i]*( (rR[i]+a(times)) * (1 - R[i]/Ki[i]) - (t(cC[i,])%*%C))
					#dR[i] = R[i]*( (rR[i]) * (1 - R[i]/Ki[i]) - (t(cC[i,])%*%C))

				}

				###Consumer dynamics: LV consumption
				dC = C 
				for( i in 1:nCsp){
					dC[i] = C[i] * ( rC[i] *(eFc[i]*cC[,i])%*%R -(t(cP[i,])%*%P)- muC[i] )
				}

				###Predator dynamics: LV consumption
				dP = P 
				for( i in 1:nPsp){
					dP[i] = P[i] * ( rP[i] *(eFp[i]*cP[,i])%*%C - muP[i] )
				}

				a = a(times) 

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
					dR[i] = R[i]*( (rR[i]) * (1 - R[i]/Ki[i]) - (t(cC[i,])%*%C))
				}

				###Consumer dynamics: LV consumption
				dC = C 
				for( i in 1:nCsp){
					dC[i] = C[i] * ( rC[i] *(eFc[i]*cC[,i])%*%R -(t(cP[i,])%*%P)- muC[i] )
				}

				###Predator dynamics: LV consumption
				dP = P 
				for( i in 1:nPsp){
					dP[i] = P[i] * ( rP[i] *(eFp[i]*cP[,i])%*%C - muP[i] )
				}

			return(list(c(dR,dC,dP)))

			}
	}


	#=============================================================================
	# Run the population models.
	#=============================================================================
	winit = c(matrix(3,nspp,1))
	out=NULL
	out$out=ode(y=winit,times=times,func=food_web,parms=parms)
	out$spp_prms=parms
	return(out)

}



#=============================================================================
# Information theoretic assessment of the foodweb.
#=============================================================================
#=============================================================================
# This code adapts the approach of  Rutledge, Basore, and Mulholland 1976
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

}
#=============================================================================
# get_mI
# Average mutual information in the Rutldege web model. The MI is interpreted 
# as the uncertainty resolved by knowing food web structure.
# When energy flow is based on pop size, then the "fij" are the per food source
# conversion rate to new population (I think this should be scaled by the freq of  
# the population i.e. e*c*pop_freq vs. e*c only). 
#
# freq 			A frequency distribution that is standardized to 1
# fij 			The transition rates from the Rutledge web 
#=============================================================================

get_mI = function (freq=freq, fij=fij) {
	ncol1 = length(freq)

	mI = NULL
	mI$mean = 0
	mI$per = matrix(0,dim(fij)[1],dim(fij)[2]) #How much does each link contribute? 
	for (k in 1:ncol1){ 
		for(j in 1:ncol1){
			mI$per[k,j] = freq[k]*fij[j,k] *log(fij[j,k]/freq[j] )
			if(is.na(mI$per[k,j])){mI$per[k,j]=0 } 
		}
	}
	mI$mean = sum(rowSums(mI$per))
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
#=============================================================================

rutledge_web = function (spp_list = spp_list, pop_ts=pop_ts, spp_prms=spp_prms) {

	tuse = dim(pop_ts)[1]
	nRsp=spp_list[1] #Resource species
	nCsp=spp_list[2] #Consumer species
	nPsp=spp_list[3] #Predator species
	nspp = nRsp+nCsp+nPsp #Total number of species

	ncol1 = nRsp+2*nCsp+2*nPsp+1
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

	#Information theoretic quantities 
	rweb$sD = matrix(0,tuse, 1) #Shannon entropy
	rweb$mI_mean = matrix(0,tuse, 1) #Mutual Information
	rweb$mI_per =array(c(matrix(0,ncol1,nrow1),matrix(0,ncol1,nrow1)), dim = c(ncol1,nrow1,(tuse))) 
	rweb$ce = matrix(0,tuse, 1) #Conditional Entropy

	###For now, loop through each time step. Is there a faster way to do this with matrix math? 
	for(n in 1:(tuse)) { 

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
		diag(rweb$fijQi[1:nRsp,1:nRsp,n]) = (rR)*R1[,1]
		#Consumer production: 
		rweb$fijQi[(1+nRsp):(nRsp+nCsp),1:nRsp,n] = t((frC*R1*cC)*t(C1r) )
		#Predator production: 
		rweb$fijQi[(1+nRsp+nCsp):nspp,(1+nRsp):(nRsp+nCsp),n] = (frP*C1p*cP)*t(P1)

		##Energetic loss: 
		#Resource to consumer: 
		rweb$fijQi[(nspp+1):(nspp+nCsp),1:nRsp,n ] = t(( R1*cC)*t(C1r) ) - t((frC*R1*cC)*t(C1r) )
		#Consumer to Predator: 
		rweb$fijQi[(nspp+nCsp+1):(nspp+nCsp+nPsp),(1+nRsp):(nRsp+nCsp),n ] = t( (C1p*cP)*t(P1) - (frP*C1p*cP)*t(P1) )

		##Mortality loss:
		#Resource: 
		rweb$fijQi[ncol1,1:nRsp,n ] =(rR)*R1[,1]*(R1[,1]/Ki)
		#Consumer:
		rweb$fijQi[ncol1,(1+nRsp):(nRsp+nCsp),n ] = t((t(C1r)*fmuC)[1,])
		#Predator
		rweb$fijQi[ncol1,(1+nRsp+nCsp):nspp,n ] = (t(P1)*fmuP)[1,]

		#Now make Qi/Pi
		rweb$Qi[,n] = rowSums( rweb$fijQi[,,n]) #*as.numeric(lower.tri(fijQi[,,n])) )
		
		#Now make fij by standardizing (dividing by Qi)
		rweb$fij[,,n] = rweb$fijQi[,,n]/ matrix(rweb$Qi[,n],ncol1,nrow1,byrow=T)
		rweb$fij[,,n][is.na(rweb$fij[,,n])] = 0

		#Pi should be the same as Qi[,n+1]. Should also work equally well before or after 
		#standardizing Qi.
		#Pi[,n] = fij[,,n]%*%Qi[,n] 
		
		#Standardize Qi to be proportions: 
		rweb$Qi[,n] = rweb$Qi[,n]/sum(rweb$Qi[,n])

		#Information theoretic quantities 
		rweb$sD[n] = shannon_D(rweb$Qi[,n])
		mI_temp = get_mI (rweb$Qi[,n], rweb$fij[,,n]  )
		rweb$mI_mean[n] = mI_temp$mean
		rweb$mI_per[,,n] = mI_temp$per
		rweb$ce[n] = rweb$sD[n] - rweb$mI_mean[n]
	}
return (rweb)
}
