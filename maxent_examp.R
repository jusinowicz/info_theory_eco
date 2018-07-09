#=============================================================================
#Recovering the exponential distribution with the Max Ent approach
#Learning the Optimization library CVXR
#Try to reproduce the Maximum Entropy example from 
#http://www.di.fc.ul.pt/~jpn/r/maxent/maxent.html
#based on the cvxr Direct Standardization example from 
#https://cvxr.rbind.io/post/examples/cvxr_direct-standardization/
#=============================================================================

library(CVXR)

n = 100
a = seq(0,10,len=n) # theoretical support [0,+oo) but we assume a light tail 
lambda = .5
A = as.data.frame(a)

##Variable we're interested in
m = nrow(A) 

## Given population mean of features
b = as.matrix(1/lambda)           # f_1(n) = n, ie, E[f_1(n)] = E[n] = 1/lambda
rownames(b)="a"

## Construct the direct standardization problem
w = Variable(m)
objective = sum(entr(w)) #Max Ent constraint
constraints = list(w>=0, sum(w) == 1, t(A) %*% w == b)
prob = Problem(Maximize(objective), constraints)

## Solve for the distribution weights
result = solve(prob)
weights = result$getValue(w)


##Plotting
plot(a,weights, pch=20, ylab="")
diff1 = dexp(0,rate=lambda) /weights[1]
curve(dexp(x,rate=lambda)/diff1, col="red", add=T)

# ##Which reproduces the following: 
# library(CVXfromR)

# n = 100
# a = seq(0,10,len=n) # theoretical support [0,+oo) but we assume a light tail 
# lambda = .5

# A = matrix(a, ncol=n)
# b = 1/lambda           # f_1(n) = n, ie, E[f_1(n)] = E[n] = 1/lambda

# # ref: web.cvxr.com/cvx/examples/cvxbook/Ch07_statistical_estim/html/maxent.html
# # entr(x)=-x*log(x), elementwise entropy function [cvxr.com/cvx/doc/funcref.html]
# cvxcode = "
#     variables pmaxent(n)
#     maximize( sum(entr(pmaxent)) )
#     sum(pmaxent) == 1;
#     A * pmaxent == b;
# "

# # it takes sometime to run a matlab session
# opt.vals = CallCVX(cvxcode, const.vars=list(n=n, A=A, b=b),
#                     opt.var.names="pmaxent", 
#                     setup.dir="C:\\Users\\jpn.INFORMATICA\\Software\\_Langs\\cvx")

# plot(a,opt.vals$pmaxent, pch=20, ylab="")
# diff = dexp(0,rate=lambda) / opt.vals$pmaxent[1] # scale back to maxent approx
# curve(dexp(x,rate=lambda)/diff, col="red", add=T)