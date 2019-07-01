#=============================================================================
# Birth-death process for a single species based on the logistic equation. 
# This is mainly to check against the analytical (mean-field) solution
# The key with specifying this model as a birth-death process is that the 
# combination N*(b(N) -d(N)) = logistic. This means there are numerous equivalent 
# ways to specify the birth and death rates. Here, I choose:
#
# n(t+1) = n*(1+r*(1-n/K)),
# b(n) = b1-b2*n
# d(n) = d1+d2*n
# 
# Where b1 = r, b2 =0 , d1 = 0, d2 = r/K
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


#Birth rate
r = 1.1

#Carrying capacity
K=5

#Mean death rate
ds = r/K

#Mean immigration rate
gs = 0

#Matrix extent
Ki = 50

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

#Lower diagonal is (density dependent) death rates
diag(Q[-1,]) = ds*(diag(Qp[-1,]))^2 

#Upper diagonal is birth 
diag(Q[,-1])[2:Ki] =gs+(diag(Qp[-1,])*r)[1:(Ki-1)]


diag(Q) = -rowSums(Q)

#This is to make the matrix recurrent, in order to calculate the 
#quasi-stationary distribution (see e.g. White 1969)
Q[Ki,1] = Q[Ki,Ki+1]
Q[1,2] = r
Q[1,1] = -(r) #+ Q[Ki,Ki]

#Take the subset, which is equal to the matrix we want
Q= Q[1:Ki,1:Ki]

#=============================================================================
# Define CVXR variables that we're interested in
#=============================================================================
A = as.data.frame(1:(Ki))
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
pop_bdi_d = get_bdi_lgs2(r,ds,gs,K)

#The more traditional M/M/1 approach where the time to the next
#event (birth/death) is exponentially distributed

pop_mm1_d = get_mm1_lgs2(r,ds,gs,K)

#=============================================================================
#Check how well distributions match
#=============================================================================
q1=hist(pop_mm1_d,probability=TRUE,breaks=0:(Ki+1))
bd1 = hist(pop_bdi_d,probability=TRUE,breaks=0:(Ki+1))
plot(0:(length(q1$density)-1), q1$density/(sum(q1$density)),ylim=c(0,1))
points(0:(length(bd1$density)-1),bd1$density/(sum(bd1$density)),col="blue")
points(0:(length(weights)-1),weights,col="red")

#How well does the average population size match to the carrying capacity? 
K_act = K
mean(pop_mm1_d)
mean(pop_bdi_d)
sum(weights*1:Ki)
