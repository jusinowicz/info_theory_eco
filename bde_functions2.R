#=============================================================================
# Functions corresponding to birth-death process code 
# Includes population processes, complex birth/death processes, and functions
# related to various information theoretic quantities (maximum entropy, 
# entropy distance, entropy deficiency)
#=============================================================================



#=============================================================================
# FUNCTION DEFINITIONS
#=============================================================================

#=============================================================================
# get_bdi: Simulation of population process as birth-death-immigration process
# bs		arrival rate (birth)
# ds			service rate (death)
# nsims=30000 	number of simulations
# ip=1 			initial population 
#=============================================================================

get_bdi = function (bs, ds, nsims=30000, ip=1) {

	pop = matrix(0,nsims,1) #Population matrix (queue)
	pop[1] = ip
	p = runif(nsims,0,1) #Probabilities of birth/death
	pt =1 

	for (t in 2:nsims) { 
 		if(pt>0) {
 			pop[t] = pt
			pt=ifelse (p[t]<bs/(bs+ds),
	            pt+1, # arrival
	            pt-1) # departure
		} else {
			pop[t] = pt
			if(p[t]<bs/(bs+ds)){pt=1}
	           
		}
	}
	return(pop)
}

#=============================================================================
# get_mm1: Simulation of population process as classic M/M/1 queue 
# credit to Tanya Leise, code adapted from
# https://tleise.people.amherst.edu/Math365Spring2014/Rscripts/QueueSimulator.R
#
# l1			arrival rate (birth)
# mu1			service rate (death)
# nsims=30000 	number of simulations
# ip=1 			initial population 
#=============================================================================


get_mm1 = function (l1, mu1, nsims=30000, ip=1) {

	t = 0 # current time
	queue = 0 # start with empty queue

	# first arrival
	T1 = rexp(1,rate=l1)
	currentqueue = 1
	eventsTime = T1
	t = T1
	nEvents = 1 # total number of events that have occurred

	while (t<nsims) {
	  nEvents = nEvents+1
	  if(currentqueue>0) {
	    T1[nEvents] =rexp(1,rate=l1+mu1) # time until next event
	    # is event an arrival or departure?
	    p[nEvents] = runif(1,0,1) 
	    queue[nEvents]=currentqueue # how many have been in queue before this new event
	    currentqueue=ifelse(p[nEvents]<l1/(l1+mu1),
	                         currentqueue+1, # arrival
	                         currentqueue-1) # departure
	  } else {
	    T1[nEvents] =rexp(1,rate=l1)
	    queue[nEvents] = currentqueue
	    currentqueue = 1
	  }
	  t = t+T1[nEvents] # time of next arrival
	  eventsTime[nEvents] = T1[nEvents] # inter-event time
	}

	return(queue)

}


#=============================================================================
# get_bdi_d: Simulation of population process as birth-death-immigration process
#			 With density-dependent birth/death
# bs, bs_d		arrival rate (birth), density dependent.
# ds, ds_d		service rate (death), density dependent
# nsims=30000 	number of simulations
# ip=1 			initial population 
#=============================================================================

get_bdi_d = function (bs, ds, Ki, nsims=30000, ip=1) {

	pop = matrix(0,nsims,1) #Population matrix (queue)
	pop[1] = ip
	p = runif(nsims,0,1) #Probabilities of birth/death
	pt =1 

	for (t in 2:nsims) { 
 		if(pt>0) {
 			pop[t] = pt
 			#Make bs ands ds density-dependent
 			bs_d = bs*pop[t]*(1-pop[t]/Ki)
			ds_d = ds*pop[t]
			pt=ifelse (p[t]<bs_d/(bs_d+ds_d),
	            pt+1, # arrival
	            pt-1) # departure
		} else {
			pop[t] = pt
			if(p[t]<bs/(bs+ds)){pt=1}
	           
		}
	}
	return(pop)
}

#=============================================================================
# get_mm1_d: Simulation of population process as classic M/M/1 queue 
# credit to Tanya Leise, code adapted from
# https://tleise.people.amherst.edu/Math365Spring2014/Rscripts/QueueSimulator.R
#			 With density-dependent birth/death
# l1, l1_d		arrival rate (birth), density dependent.
# mu1, mu1_d	service rate (death), density dependent
# nsims=30000 	number of simulations
# ip=1 			initial population 
#=============================================================================


get_mm1_d = function (l1, mu1, Ki, nsims=30000, ip=1) {

	t = 0 # current time
	queue = 0 # start with empty queue

	# first arrival
	T1 = rexp(1,rate=l1)
	currentqueue = 1
	eventsTime = T1
	t = T1
	nEvents = 1 # total number of events that have occurred

	l1_d = l1
	mu1_d = mu1

	while (t<nsims) {
	  nEvents = nEvents+1
	  if(currentqueue>0) {
	    T1[nEvents] =rexp(1,rate=l1_d+mu1_d) # time until next event
	    # is event an arrival or departure?
	    p[nEvents] = runif(1,0,1) 
	    queue[nEvents]=currentqueue # how many have been in queue before this new event
	    l1_d = l1*queue[nEvents]*(1-queue[nEvents]/Ki)
		mu1_d = mu1*queue[nEvents]
	    currentqueue=ifelse(p[nEvents]<l1_d/(l1_d+mu1_d),
	                         currentqueue+1, # arrival
	                         currentqueue-1) # departure
	  } else {
	    T1[nEvents] =rexp(1,rate=l1)
	    queue[nEvents] = currentqueue
	    currentqueue = 1
	  }
	  t = t+T1[nEvents] # time of next arrival
	  eventsTime[nEvents] = T1[nEvents] # inter-event time
	}

	return(queue)

}