#=============================================================================
# Birth-death process for a single species based on the M/M/1 qeueing process.
# Biologically, this is a birth-death-immigration process because 
# the population will tend to 0, but new individuals still appear when there is
# 0 population.
# 
# Although this model is a bit biologically weird, it is a starting place due to its
# analytical tractability and known solutions. 
#
# Stationary distribution solved numerically using maximum entropy.
#
# This reproduces a known example with intraspecific density-dependent birth
# through a carrying capacity
#=============================================================================
#Load libraries
#=============================================================================
library(fitdistrplus)
library(CVXR)
source("./birth_deathR/bde_functions_lgs.R")
#=============================================================================

#Number of generations:
n=100
N= 0:n


#Mean birth rate
bs = 30

#Mean death rate
ds = 15

#Mean immigration rate
gs = 0

#Carrying capacity
K=15

#Matrix extent (just needs to be bigger than the actual carrying capacity)
Ki=80

#=============================================================================
# Using maximum entropy (through CVXR) to solve for stationary distribution
# This uses the constraints: 
#		sum(pr) = 1   (probabilities standardize to 1)
#		(t(Q)%*%pr) == 0matrix (Each transition = 0, steady-state dynamics)
#=============================================================================

#=============================================================================
#Build Q, the transition rate matrix. Each entry is the birth/death probability
# multiplied by the number of individuals. 
#=============================================================================

#Make the new extent equal to Ki
Qp = matrix( (0:(Ki+1)), (Ki+2),(Ki+2))	
Qp2 = matrix( 1:(Ki+2), (Ki+2),(Ki+2))
Qp[upper.tri(Qp)] = Qp2[upper.tri(Qp2)] 
Q = matrix(0,(Ki+2),(Ki+2))

#Upper diagonal is (density dependent) birth + immigration
diag(Q[,-1])[2:(Ki+1)] =gs+(diag(Qp[-1,])*bs*(1- (diag(Qp[-1,]))/K))[1:(Ki)]

#Lower diagonal is death rates
#Note: the density from the lower diagonal is used just to match the final
#vector size. Could also use diag(Qp)[2:Ki]
diag(Q[-1,]) =diag(Qp[-1,])*ds

diag(Q) = -rowSums(Q)

#This is to make the matrix recurrent, in order to calculate the 
#quasi-stationary distribution (see e.g. White 1969)
Q[Ki,1] = Q[Ki,Ki+1]
Q[1,2] = bs
Q[1,1] = -(bs) #+ Q[Ki,Ki]

#Take the subset, which is equal to the matrix we want
Q= Q[1:Ki,1:Ki]

#=============================================================================
# Define CVXR variables that we're interested in
#=============================================================================
A = as.data.frame(1:Ki)
m = nrow(A)

# Construct the direct standardization problem
w = Variable(m)
objective = sum(entr(w)) #Max Ent constraint

constraints = list(sum(w) == 1, ( t(Q)%*%(w)) == matrix(0,(Ki),1) )
#constraints = list(w>=0, sum(w) == 1, t(A) %*% w == b) #Classic mean constraint
prob = Problem(Maximize(objective), constraints)

## Solve for the distribution weights
#result = solve(prob,solver="SCS") #You may need to change the solver
result = solve(prob) #You may need to change the solver

weights = result$getValue(w)

#=============================================================================
#Population process, as a check
#=============================================================================
pop_bdi_d = get_bdi_d(bs,ds,gs,K)

#The more traditional M/M/1 approach where the time to the next
#event (birth/death) is exponentially distributed

pop_mm1_d = get_mm1_d(bs,ds,gs,K)

#=============================================================================
#Check how well distributions match
#=============================================================================

q1=hist(pop_mm1_d,probability=TRUE,breaks=0:(Ki+1))
bd1 = hist(pop_bdi_d,probability=TRUE,breaks=0:(Ki+1))
plot(0:(length(q1$density)-1), q1$density/(sum(q1$density)),ylim=c(0,1))
points(0:(length(bd1$density)-1),bd1$density/(sum(bd1$density)),col="blue")
points(0:(length(weights)-1),weights,col="red")

#How well does the average population size match to the carrying capacity? 
K_act = (bs-ds)*K/bs
mean(pop_mm1_d)
mean(pop_bdi_d)

#Find the mean of the probability distribution. If thesea are not equal, thne 
#something odd is going on (Most likely the matrix extent)
#sum(weights*1:Ki)
sum(weights*0:(Ki-1))
which(weights==max(weights))-1


