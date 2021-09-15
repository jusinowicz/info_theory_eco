This collection of files was used to test ideas for and produce figures in “The fitness value of ecological information in a variable world” by J. Usinowicz and M.I. O'Connor, In Prep. The files herein are focused on three discrete-time and one continuous time model that range in order of complexity. Conceptually, each model tracks and analyzes population dynamics in a fluctuating environment. 

Common code: 
Common code is used to produce random environments with a user-specified underlying distribution of states, and produce species-specific preferences for environment types by mapping a Gaussian distribution of fitness onto the random environment distribution. See "env_functions.R" for details. 

Common code is used to provide information theoretic (IT) analyses of environments and population dynamics. The core IT metrics used in these analyses are the Shannon entropy, the conditional information, and the mutual information. See "info_theory_functions.R"

Population models: 
1. The simplest model in "simple_fvoi.R" is a simple multiplicative growth process where all resoruces are "spent" each generation and there is no survival. It corresponds with the classic Kelly betting problem/optimum. See e.g. Cover and Thomas, 2006, Elements of Information Theory for a great overview of this model and its significance in probability theory/IT. This model also corresponds to the model employed in Donaldson-Matasci et al. 2010.  

This code reproduces the model described in Box 1 in the main manuscript and is used to produce Figure 2A and 2C.  

2. The next model in "simple_dorm_fvoi.R" adds in the possibility to retain resources each generation, which also means that a portion of the population can survive each generation. In probability betting problems this is the idea of bet-hedging. Biologically, this maps onto the concept of dormancy and has been used by a number of authors to investigate the dynamics of annual plants and freshwater crustaceans (e.g. Daphnia and copepods) that maintain seed/egg banks. 

This code reproduces the model described in Box 1 in the main manuscript and is used to produce Figure 2B and 2D. 

3. The model in "lott_info_inv.R" is the lottery model (e.g. Chesson and Warner, 1981), modified here to include the possibility for emergence to use information. In particular, we allow the environment to produce a cue, and populations use this cue to decide what proportion of individuals to emerge and reproduce. This is done by generating conditional probability tables between cue and environment states. In terms of population dynamics, the lottery model is similar to the dormancy model but introduces density-dependent competition. All emerging individuals capture resources in proportion to their numbers each generation. 

This code reproduces the model described in Box 2 in the main manuscript and is used to produce Figure 3. 

4. The code in "lott_info_inv_allcomp.R" is the workhorse implementation of the lottery model with information used to generate many scenarios of interspecific competition and keep track of all of the analyses for each scenario. Specifically, we used this code to generate many iterations across different scenarios of resource use overlap between two competing species. Increasing levels of resource overlap increase competition and reduce the potential for coexistence. 

This is the code used to produce the data in Figure 4. 

5. The code in "cc_fvoi.R" is the Lotka-Volterra model of competition from Gil et al. 2018. It includes some code to double check that the analytical solution to the FVOI for this model works. 

This is the code used to produce Figure 5. 

Figures: 
The code in "figures.R" is used to make all of the figures in the manuscript. It loads stored data files, all of which are save with a ".var" ending. This is the map between files and figures: 

load("./data/fvoi_plot1.var") #Figure 1,3,5
#"ni_simple.var" #Figure 2
#"env_fit2.var"# #Figure 4
# "dm_simp.var" #Figure 2 


