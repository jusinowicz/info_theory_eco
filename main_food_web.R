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
# Define parameters 
#=============================================================================

###Time to simulate over: 
tend = 1000 #Length of numerical simulation
delta1 = 0.01 #Time increments
times  = seq(from = 0, to = tend, by = delta1)
tl = length(times)

###Number of species at each trophic level: 
nRsp=2 #Resource species
nCsp=3 #Consumer species
nPsp=1 #Predator species
nspp = nRsp+nCsp+nPsp #Total number of species

###Parameters of the model:
#Resource: Nearly identical resource dynamics: 
rR = matrix(rnorm(nRsp,30,0.1), nRsp, 1) #intrinsic growth
Ki = matrix(rnorm(nRsp,50,0.1), nRsp, 1) #carrying capacity

#Random resources:
c = 0.1
#Make the a variable -- see the documentation for forcings for an example
amp = 1
xint = 0
a = approxfun( x = times, y = amp*exp(rnorm(times)+xint), method = "linear", rule = 2) 
a_t = exp(amp*rnorm(times)+xint)
#a = approxfun( x = times, y = amp*rnorm(times)+xint, method = "linear", rule = 2) 
#a_t = amp*rnorm(times)+xint
print( paste("Var in a(t) = ", var(a_t),sep="")) 
a_m = mean(a_t)
vara = var(a_t)


#Consumers: 
rC = matrix(rnorm(nCsp,.8,0.1), nCsp, 1) #intrisic growth
eFc = matrix(1,nCsp,nRsp) # just make the efficiency for everything 1 for now
muC = matrix(rnorm(nCsp,0.6,0.1), nCsp, 1) #mortality rates
#Consumption rates: 
#This is one way to generate a hierarchy where each species 
#predominantly feeds on particular resource. 
dspp = abs((nCsp - nRsp))
hier1 = c(1, 1/2, 1/3) #Sum of this is total consumption rate when efC=1
hier2 = c(1/3,1/3,1/3)
cC = hier2 
for( n in 1:(nCsp-dspp)) {
	cC = cbind(cC, shifter(hier1,n))
}
cC = as.matrix(cC[1:nRsp,1:nCsp ])

#Predators: 
rP = matrix(rnorm(nPsp,0.8,0.1), nPsp, 1) #intrisic growth
eFp = matrix(1,nPsp,nCsp) # just make the efficiency for everything 1 for now
muP = matrix(rnorm(nPsp,0.6,0.1), nPsp, 1) #mortality rates
#Consumption rates: 
#This is one way to generate a hierarchy where each species 
#predominantly feeds on particular resource. 
dspp = ((nPsp - nCsp))
if(dspp<0){dspp = 0 }
hier1 = matrix(c(3/4, 1/3, 1/4, 1/6),4,1) #Sum of this is total consumption rate when efC=1
cP = hier1
for( n in 1:(nPsp-dspp-1)) {
	cP = cbind(cP, shifter(hier1,n))
}
cP = as.matrix(cP[1:nCsp,1:nPsp])


#Pass all of these parameters as a list
parms = list(nspp=nspp, nRsp = nRsp, nCsp = nCsp, nPsp =nPsp,
	rR = rR, Ki =Ki,
	rC = rC, eFc = eFc, muC = muC, cC = cC,
	rP = rP, eFp = eFp, muP = muP, cP = cP,
	a = a
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

	#a = a(times) 

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
plot(out[,"1"],t="l",col="red",ylim = c(0,max(out[tl,2:nspp],na.rm=T)))
for( n in 2:(nRsp) ) {
lines(out[,paste(n)],t="l",col="red")
}
for( n in ( (nRsp+1):(nRsp+nCsp) ) ) {
lines(out[,paste(n)],t="l",col="blue")
}
for( n in ((nRsp+nCsp):nspp ) ) {
lines(out[,paste(n)],t="l")
}

#=============================================================================
# Information theoretic assessment of the foodweb.
#=============================================================================
#=============================================================================
# This section is as per Rutledge, Basore, and Mulholland 1976
#=============================================================================
## This code takes the ODEs and converts them to a biomass balance matrix and 
## transition matrix. 
## This version creates a compartment for each "event" where biomass is gained 
## or loss. This includes birth, death, and "inefficiency" in the form of the 
## way that biomass consumed translates to new population biomass. 

ts = tl
tbegin = tl/2
pop_ts = out[tbegin:(ts),2:(nspp+1)]
tuse = dim(pop_ts)[1]
ncol1 = nRsp+2*nCsp+2*nPsp+1
nrow1 = ncol1

###Generate the quantities that describe "energy" (biomass?) flow through
###food web.
fij = array(c(matrix(0,ncol1,nrow1),matrix(0,ncol1,nrow1)), dim = c(ncol1,nrow1,(tuse))) 
Qi = matrix(0,ncol1,(tuse))
fijQi = array(c(matrix(0,ncol1,nrow1),matrix(0,ncol1,nrow1)), dim = c(ncol1,nrow1,(tuse))) 
#Pi = matrix(0,ncol1,(tuse))

###For now, loop through each time step. Is there a faster way to do this with matrix math? 
for(n in 1:(tuse)) { 

	#Make a block matrix of biomass transfer for each trophic level 

	#Timestep 1
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
	diag(fijQi[1:nRsp,1:nRsp,n]) = (rR+a_t[(n)])*R1[,1]
	#Consumer production: 
	fijQi[(1+nRsp):(nRsp+nCsp),1:nRsp,n] = t((frC*R1*cC)*t(C1r) )
	#Predator production: 
	fijQi[(1+nRsp+nCsp):nspp,(1+nRsp):(nRsp+nCsp),n] = (frP*C1p*cP)*t(P1)

	##Energetic loss: 
	#Resource to consumer: 
	fijQi[(nspp+1):(nspp+nCsp),1:nRsp,n ] = t(( R1*cC)*t(C1r) ) - t((frC*R1*cC)*t(C1r) )
	#Consumer to Predator: 
	fijQi[(nspp+nCsp+1):(nspp+nCsp+nPsp),(1+nRsp):(nRsp+nCsp),n ] = t( (C1p*cP)*t(P1) - (frP*C1p*cP)*t(P1) )

	##Mortality loss:
	#Resource: 
	fijQi[ncol1,1:nRsp,n ] =(rR+a_t[(n)])*R1[,1]*(R1[,1]/Ki)
	#Consumer:
	fijQi[ncol1,(1+nRsp):(nRsp+nCsp),n ] = t((t(C1r)*fmuC)[1,])
	#Predator
	fijQi[ncol1,(1+nRsp+nCsp):nspp,n ] = (t(P1)*fmuP)[1,]

	#Now make Qi/Pi
	Qi[,n] = rowSums( fijQi[,,n]) #*as.numeric(lower.tri(fijQi[,,n])) )
	
	#Now make fij by standardizing (dividing by Qi)
	fij[,,n] = fijQi[,,n]/ matrix(Qi[,n],ncol1,nrow1,byrow=T)
	fij[,,n][is.na(fij[,,n])] = 0

	#Pi should be the same as Qi[,n+1]. Should also work equally well before or after 
	#standardizing Qi.
	#Pi[,n] = fij[,,n]%*%Qi[,n] 
	
	#Standardize Qi to be proportions: 
	Qi[,n] = Qi[,n]/sum(Qi[,n])

}

#Generate quantities for the maximum entropy distribution, i.e. uniform: 
pop_me = runif(nspp)
me_freq = pop_me/matrix(sum(pop_me),length(pop_me),1)

###"Thoroughput " diversity (distribution of energetic flow through web)
D_pop1 = - colSums ( Qi*log(Qi),na.rm=T )
D_me = - sum ( me_freq*log(me_freq) )

###Average mutual information -- uncertainty resolved by knowing food web structure
#When energy flow is based on pop size, then the "fij" are the per food source
#conversion rate to new population (I think this should be scaled by the freq of  
#the population i.e. e*c*pop_freq vs. e*c only). 
mI = matrix (0,tuse,1)
for(n in 2:(tuse)) { 
	for (k in 1:ncol1){ 
		for(j in 1:ncol1){

			ntemp = Qi[k,(n-1)]*fij[j,k,(n-1)] *log(fij[j,k,(n-1)]/Qi[j,n] )
			if(is.na(ntemp)){ntemp =0 }
			mI[(n-1)]= mI[( n-1)]+ntemp
		}

	}
}

###Conditional entropy -- remaining uncertainty
S_ce = D_pop1 - mI

#=============================================================================
#### Compare to this version of the matrix, where the compartments are just per-
#### species

fij_s = array(c(matrix(0,ncol1,nrow1),matrix(0,ncol1,nrow1)), dim = c(ncol1,nrow1,(tuse))) 
Qi_s = matrix(0,ncol1,(tuse))
fijQi_s = array(c(matrix(0,ncol1,nrow1),matrix(0,ncol1,nrow1)), dim = c(ncol1,nrow1,(tuse))) 

for(n in 1:(tuse)) {

#fijQi_s = 

}