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
bs = 1.1

#Mean death rate
ds = 1

#Mean immigration rate
gs = 0

#Carrying capacity
Ki=5

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

#Make the new extent equal to carrying capacity, Ki
Qp = matrix( (0:Ki), (Ki+1),(Ki+1))
Qp2 = matrix( 1:(Ki+1), (Ki+1),(Ki+1))
Qp[upper.tri(Qp)] = Qp2[upper.tri(Qp2)] 
Q = matrix(0,(Ki+1),(Ki+1))

#Fix the first 0 entries for immigration
Q[1,2] = gs
Q[1,1] = -gs

#Lower diagonal is death rates
#Note: the density from the lower diagonal is used just to match the final
#vector size. Could also use diag(Qp)[2:Ki]
diag(Q[-1,]) =diag(Qp[-1,])*ds

#Upper diagonal is (density dependent) birth + immigration
diag(Q[,-1])[2:Ki] =gs+(diag(Qp[-1,])*bs*(1- (diag(Qp[-1,]))/Ki))[1:(Ki-1)]

#Main diagonal is the -(bs+ds) because system will not stay in current state
diag(Q)[2:(Ki+1)] = -( diag(Q[-1,])+c(diag(Q[,-1])[2:Ki],0))



#=============================================================================
# Define CVXR variables that we're interested in
#=============================================================================
A = as.data.frame(0:Ki)
m = nrow(A)

# Construct the direct standardization problem
w = Variable(m)
objective = sum(entr(w)) #Max Ent constraint

constraints = list(sum(w) == 1, ( t(Q)%*%(w)) == matrix(0,(Ki+1),1) )
#constraints = list(w>=0, sum(w) == 1, t(A) %*% w == b) #Classic mean constraint
prob = Problem(Maximize(objective), constraints)

## Solve for the distribution weights
result = solve(prob,solver="SCS") #You may need to change the solver
weights = result$getValue(w)

#=============================================================================
#Population process, as a check
#=============================================================================
pop_bdi_d = get_bdi_d(bs,ds,gs,Ki)

#The more traditional M/M/1 approach where the time to the next
#event (birth/death) is exponentially distributed

pop_mm1_d = get_mm1_d(bs,ds,gs,Ki)

#=============================================================================
#Check how well distributions match
#=============================================================================

q1=hist(pop_mm1_d,probability=TRUE,breaks=0:(Ki+1))
bd1 = hist(pop_bdi_d,probability=TRUE,breaks=0:(Ki+1))
plot(0:(length(q1$density)-1), q1$density/(sum(q1$density)),ylim=c(0,1))
points(0:(length(bd1$density)-1),bd1$density/(sum(bd1$density)),col="blue")
points(0:(length(weights)-1),weights,col="red")
