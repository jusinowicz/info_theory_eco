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
# Define parameters 
#=============================================================================

###Time to simulate over: 
tend = 1000 #Length of numerical simulation
times  = seq(from = 0, to = tend, by = 0.01)

###Number of species at each trophic level: 
nRsp=3 #Resource species
nCsp=4 #Consumer species
nPsp=3 #Predator species
nspp = nRsp+nCsp+nPsp #Total number of species

###Parameters of the model:
#Resource: Nearly identical resource dynamics: 
rR = matrix(rnorm(nRsp,20,0.1), nRsp, 1) #intrinsic growth
Ki = matrix(rnorm(nRsp,50,0.1), nRsp, 1) #carrying capacity

#Consumers: 
rC = matrix(rnorm(nCsp,3,0.1), nCsp, 1) #intrisic growth
eFc = matrix(1,nCsp,1) # just make the efficiency for everything 1 for now
muC = matrix(rnorm(nCsp,0.6,0.1), nCsp, 1) #mortality rates
#Consumption rates: 
#This is one way to generate a hierarchy where each species 
#predominantly feeds on particular resource. 
dspp = abs((nCsp - nRsp))
hier1 = c(1, 1/2, 1/3) #Sum of this is total consumption rate when efC=1
hier2 = c(1/3,1/3,1/3)
cC = hier2 
for( n in 1:(nCsp-dspp)) {
	cC = rbind(cC, shifter(hier1,n))
}

#Predators: 
rP = matrix(rnorm(nPsp,1,0.1), nPsp, 1) #intrisic growth
eFp = matrix(1,nPsp,1) # just make the efficiency for everything 1 for now
muP = matrix(rnorm(nPsp,0.6,0.1), nPsp, 1) #mortality rates
#Consumption rates: 
#This is one way to generate a hierarchy where each species 
#predominantly feeds on particular resource. 
dspp = ((nPsp - nCsp))
if(dspp<0){dspp = 0 }
hier1 = c(3/4, 1/3, 1/4, 1/6) #Sum of this is total consumption rate when efC=1
cP = hier1
for( n in 1:(nPsp-dspp-1)) {
	cP = rbind(cP, shifter(hier1,n))
}

#Pass all of these parameters as a list
parms = list(nspp=nspp, nRsp = nRsp, nCsp = nCsp, nPsp =nPsp,
	rR = rR, Ki =Ki,
	rC = rC, eFc = eFc, muC = muC, cC = cC,
	rP = rP, eFp = eFp, muP = muP, cP = cP
 )

#=============================================================================
# Define population dynamics
#=============================================================================
### Function specifying the full dynamics of all trophic levels:
food_web = function(times,sp,parms){
	nspp = parms$nspp 
     R = matrix(sp[1:nRsp], nRsp, 1)
     C = matrix(sp[(nRsp+1):(nRsp+nCsp)], nCsp, 1)
     P = matrix(sp[(nRsp+nCsp+1):nspp], nPsp, 1)


	###Resource dynamics: Logistic growth, reduced by consumption
	dR = R
	for( i in 1:nRsp){
		dR[i] = rR[i] *R[i] * ((1 - R[i])/Ki[i]) - (t(cC[,i])%*%C)*R[i]
	}

	###Consumer dynamics: LV consumption
	dC = C 
	for( i in 1:nCsp){
		dC[i] = rC[i] *C[i] * ( (eFc[i]*cC[i,])%*%R - muC[i] )
	}

	###Predator dynamics: LV consumption
	dP = P 
	for( i in 1:nPsp){
		dP[i] = rP[i] *P[i] * ( (eFp[i]*cP[i,])%*%C - muP[i] )
	}

return(list(c(dR,dC,dP)))

}

#=============================================================================
# Run the population models.
#=============================================================================
winit = c(matrix(3,nspp,1))
out=ode(y=winit,times=times,func=food_web,parms=parms)
#=============================================================================
# Basic plots of species populations. 
#=============================================================================
plot(out[,"1"],t="l",col="red",ylim=c(0,1))
for( n in 2:(nRsp) ) {
lines(out[,paste(n)],t="l",col="red")
}
for( n in ( (nRsp+1):(nRsp+nCsp) ) ) {
lines(out[,paste(n)],t="l",col="blue")
}
for( n in ((nRsp+nCsp):nspp ) ) {
lines(out[,paste(n)],t="l")
}




