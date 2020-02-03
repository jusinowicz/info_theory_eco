# info_theory_eco
Using information theory to solve ecology

This repository is in its infancy. The ultimate goal is to solve ecological and evolutionary problems using tools from information theory that are currently applied to a range of problems in other disciplines. The starting focus is on using maximum entropy and entropy distance approaches to understand population dynamics, and the dynamics of competing species in particular. 

The main thrust of the work represented in this repository is currently contained in the folders "random_foodwebs" and "scenario_foodwebs." These two folders are based off of the same general code and approach which includes: 

1. Generate a  food web with a generic underlying dynamic model. This model is currently based on the model found in Gilbert et al. 2014 ( Gilbert, B., T. D. Tunney, K. S. McCann, J. P. DeLong, D. A. Vasseur, V. Savage, J. B. Shurin, A. I. Dell, B. T. Barton, and C. D. Harley. 2014. A bioenergetic framework for the temperature dependence of trophic interactions. Ecology Letters 17:902–914). Food-web includes resource, herbivore, and predator: 
	A. Resource is based on a consumer-resource model, with added predators 
		1. Competition between consumers and resources emerges from consumption
		2. Parameters at each level can be made a function of temperature. 
	 B. Resources can be stochastic due to environmental fluctuations. 
	 C. Relative non-linearity allows 2 consumers per Resource
2. Either Generate a bunch of random food webs (random_) or build a particular set (scenario_)
3. Use dynamic information theoretic metrics to understand the resulting food-web structures. 
4. Visualize results from both the foodweb and its information theoretic properties 

In additon, this repository also contains:
In bd_models: 

1. Reworked examples of maximum entropy approaches for numerically solving best-fit distributions subject to (Lagrangian) constraints.
  
  *  max_ent.R: The only goal of this script is to start learning CVXR. Here, I reproduce the Maximum Entropy example from http://www.di.fc.ul.pt/~jpn/r/maxent/maxent.html based on the CVXR Direct Standardization example from https://cvxr.rbind.io/post/examples/cvxr_direct-standardization  
  
2. An example of a birth-death process where births are density dependent (based on carrying capacity). This reproduces an example from https://tice.agroparistech.fr/coursenligne/courses/MODULEINTEGRATIF17MO/
But also builds this into a more general approach to reproduce classic solutions to the deterministic logistic equation. 

  * bd_density_entropy.R: This file reproduces the birth-death process with density dependence using simulations, and by solving for the stationary distribution with a maximum entropy approach. The maximum entropy approach uses the description of equilibrium in the matrix of transition rates as a constraint and solves for it numerically. The numerical solving is done using the CVXR package. 
  * bd_density_entropy_quasi2: Almost the same as bd_density_entropy.R, but where the quasi-stationary distribution has been solved instead. The mean of the quasi-stationary distribution corresponds to the classic solution of the deterministic logistic equation. 
  * logistic_entropy2: The same as bd_density_entropy_quasi, but where a different form of the logistic model has been implemented. 
  * bde_functions_lgs.R: The code for the population simulations. There are several versions of the function, casting it  as a more traditional M/M/1 queing simulation with exponentially distributed wait times, and as a density-dependent logistic growth process. 
