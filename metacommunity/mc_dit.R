#=============================================================================
#Load libraries
#=============================================================================
library(tidyverse)
library(RandomFields)
library(vegan)

#Patrick Thompson's metacommunity model functions: 
source("./mcomsimr-master/R/MC_simulate.R")
source("./mcomsimr-master/R/sim_setup.R")
source("../info_theory_functions/info_theory_functions.R")

#=============================================================================



###This is an example of just running the model. 
# a1=simulate_MC(patches = 16, species = 1, env1Scale = 100, dispersal = 0.001, min_inter = 1, max_inter = 1, env_niche_breadth = 1)
# env.df = a1$env.df
# g<-ggplot(env.df, aes(x = time, y = env1, group = patch, color = factor(patch)))+
#       geom_line()+
#       scale_color_viridis_d(guide=F)
# print(g)



### Remove any temporal variance in the landscape (cannot figure out how to do this with RandomFields): 
patches = 16
nspp = 2
landscape = landscape_generate(patches = patches, plot = TRUE)

var_b = 1
env1Scale = 100
timesteps = 1200 
burn_in = 1200  
initialization = 200

#Make the env.df outside of the simulate_MC function with RandomFields.
#This code makes a single time realization of the RandomFields model, then just repeats it over the right number of 
#time steps (timesteps + burn_in)
model = RMexp(var=var_b, scale=env1Scale) + # with variance 4 and scale 10
       RMnugget(var=0) + # nugget
        RMtrend(mean=0.01) # and mean
RF = RFsimulate(model = model, x = landscape$x*10, y = landscape$y*10, spConform=FALSE)
env.df = data.frame(env1 = decostand(RF,method = "range"), patch = 1:nrow(landscape), time = rep(1:(timesteps+burn_in), each = nrow(landscape)))
env.initial = env.df[env.df$time == 1,]

#I've added the "constant = TRUE" option to the function to remove demographic stochasticity. 
#Added mortality to the model. With m = 0, this is the original model. This was necessary to get extinction without
#stochasticity. 
a1=simulate_MC(patches = patches , species = nspp, env1Scale = env1Scale, dispersal = .000001, min_inter = 1.1, max_inter = 1.1,
	 env_niche_breadth = 1, max_r = c(5,4), m = 0.5, env.df = env.df, landscape = landscape, constant = TRUE)

#Plot the environment by patch 
g<-ggplot(a1$dynamics.df, aes(x = time, y = env, group = patch, color = factor(patch)))+
      geom_line()+
      scale_color_viridis_d(guide=F)
print(g)
