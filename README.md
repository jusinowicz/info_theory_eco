# info_theory_eco
Using information theory to solve ecology

This repository is in its infancy. The ultimate goal is to solve ecological and evolutionary problems using tools from information theory that are currently applied to a range of problems in other disciplines. The starting focus is on using maximum entropy and entropy distance approaches to understand population dynamics, and the dynamics of competing species in particular. 

This repository currently contains:

1. Reworked examples of maximum entropy approaches for numerically solving best-fit distributions subject to (Lagrangian) constraints.
  
  *  max_ent.R: The only goal of this script is to start learning CVXR. Here, I reproduce the Maximum Entropy example from http://www.di.fc.ul.pt/~jpn/r/maxent/maxent.html based on the CVXR Direct Standardization example from https://cvxr.rbind.io/post/examples/cvxr_direct-standardization  
  
2. An example of a birth-death process where births are density dependent (based on carrying capacity). This reproduces an example from https://tice.agroparistech.fr/coursenligne/courses/MODULEINTEGRATIF17MO/
But also builds this into a more general approach to reproduce classic solutions to the deterministic logistic equation. 

  * bd_density_entropy.R: This file reproduces the birth-death process with density dependence using simulations, and by solving for the stationary distribution with a maximum entropy approach. The maximum entropy approach uses the description of equilibrium in the matrix of transition rates as a constraint and solves for it numerically. The numerical solving is done using the CVXR package. 
  * bd_density_entropy_quasi2: Almost the same as bd_density_entropy.R, but where the quasi-stationary distribution has been solved instead. The mean of the quasi-stationary distribution corresponds to the classic solution of the deterministic logistic equation. 
  * logistic_entropy2: The same as bd_density_entropy_quasi, but where a different form of the logistic model has been implemented. 
  * bde_functions_lgs.R: The code for the population simulations. There are several versions of the function, casting it  as a more traditional M/M/1 queing simulation with exponentially distributed wait times, and as a density-dependent logistic growth process. 
