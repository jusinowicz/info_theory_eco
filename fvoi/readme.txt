This collection of files was used to test ideas for and produce figures in “The fitness value of ecological information in a variable world” by J. Usinowicz and M.I. O'Connor, In Prep. The files herein are focused on three discrete-time and one continuous time model that range in order of complexity. Conceptually, each model tracks and analyzes population dynamics in a fluctuating environment. 

Common code: 
Common code is used to produce random environments with a user-specified underlying distribution of states, and produce species-specific preferences for environment types by mapping a Gaussian distribution of fitness onto the random environment distribution. See "env_functions.R" for details. 

Common code is used to provide information theoretic (IT) analyses of environments and population dynamics. The core IT metrics used in these analyses are the Shannon entropy, the conditional information, and the mutual information.

Population models: 
The simplest model in "simple_fvoi.R" is a simple multiplicative growth process, and corresponds with the classic Kelly betting problem/optimum. See e.g. Cover and Thomas, 2006, Elements of Information Theory for a great overview of this model and its significance in probability theory/IT. This model also corresponds to the model employed in Donaldson-Matasci et al. 2010. 

 