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
tl = length(times)

###Number of species at each trophic level: 
nRsp=3 #Resource species
nCsp=4 #Consumer species
nPsp=3 #Predator species
nspp = nRsp+nCsp+nPsp #Total number of species

###Parameters of the model:
#Resource: Nearly identical resource dynamics: 
rR = matrix(rnorm(nRsp,30,0.1), nRsp, 1) #intrinsic growth
Ki = matrix(rnorm(nRsp,50,0.1), nRsp, 1) #carrying capacity

#Consumers: 
rC = matrix(rnorm(nCsp,5,0.1), nCsp, 1) #intrisic growth
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
	cC = rbind(cC, shifter(hier1,n))
}

#Predators: 
rP = matrix(rnorm(nPsp,1,0.1), nPsp, 1) #intrisic growth
eFp = matrix(1,nPsp,nCsp) # just make the efficiency for everything 1 for now
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
		dR[i] = rR[i] *R[i] * (1 - R[i]/Ki[i]) - (t(cC[,i])%*%C)*R[i]
	}

	###Consumer dynamics: LV consumption
	dC = C 
	for( i in 1:nCsp){
		dC[i] = rC[i] *C[i] * ( (eFc[i]*cC[i,])%*%R -(t(cP[,i])%*%P)- muC[i] )
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
# This first section is as per Rutledge, Basore, and Mulholland 1976
#=============================================================================
# Total energy at time ts
ts = tl
pop_ts1 = out[ts,2:(nspp+1)]
pop_conv = matrix(1,nspp,1) #energetic content per individual. 
pop_tot = sum(pop_ts1*e_conv)
pop_freq1 = pop_ts1/matrix(e_tot,length(pop_ts1),1) #Energy distribution

#Generate quantities for the maximum entropy distribution, i.e. uniform: 
pop_me = runif(nspp)
me_freq = pop_me/matrix(sum(pop_me),length(pop_me),1)

###Generate the quantities that describe "energy" (biomass?) flow through
###food web.

fij = matrix(0,nspp+1,nspp+1)
fijQi = matrix(0,nspp+1,nspp+1)

#Make a block matrix of freq for each trophic level 
fR1 = matrix(pop_freq1[1:nRsp],nCsp,nRsp,byrow=T)
fC1 = matrix(pop_freq1[(nRsp+1):(nRsp+nCsp)],nPsp,nCsp,byrow=T)
fP1 = matrix(pop_freq1[(nRsp+nCsp+1):(nspp)],nCsp,nPsp,byrow=T)

frC = matrix(rC,nCsp,nRsp)
frP = matrix(rP,nPsp,nCsp)

fmuC = matrix(muC,nCsp,nRsp)
fmuP = matrix(muP,nPsp,nCsp)

#Resource
diag(fijQi[1:nRsp,1:nRsp]) = rR*(1/pop_tot - fR1[1,])/Ki*fR1[1,]
#fij[nspp+1,1:nRsp] = colSums()
#Consumer
fijQi[(1+nRsp):(nRsp+nCsp),1:nRsp] = (fR1*cC*t(fC1)) 
fijQi[nspp+1,(1+nRsp):(nRsp+nCsp) ] = (frC*t(fC1)* (fmuC/pop_tot))[,1]
#Predator
fijQi[(1+nRsp+nCsp):nspp,(1+nRsp):(nRsp+nCsp)] = (fC1*cP*t(fP1))
fijQi[nspp+1,(1+nRsp+nCsp):nspp ] = (frP*t(fP1)*fmuP/pop_tot)[,1]

#Resource
diag(fij[1:nRsp,1:nRsp]) = rR*(1/pop_tot - fR1[1,])/Ki
#fij[nspp+1,1:nRsp] = colSums()
#Consumer
fij[(1+nRsp):(nRsp+nCsp),1:nRsp] = (cC*t(fC1)) 
fij[nspp+1,(1+nRsp):(nRsp+nCsp) ] = (frC* (fmuC/pop_tot))[,1]
#Predator
fij[(1+nRsp+nCsp):nspp,(1+nRsp):(nRsp+nCsp)] = (cP*t(fP1))
fij[nspp+1,(1+nRsp+nCsp):nspp ] = (frP*fmuP/pop_tot)[,1]

###"Thoroughput " diversity (distribution of energetic flow through web)
D_pop1 = - sum ( pop_freq1*log(pop_freq1) )
D_me = - sum ( me_freq*log(me_freq) )

###Average mutual information -- uncertainty resolved by knowing food web structure
#When energy flow is based on pop size, then the "fij" are the per food source
#conversion rate to new population (I think this should be scaled by the freq of  
#the population i.e. e*c*pop_freq vs. e*c only). 


###Conditional entropy -- remaining uncertainty