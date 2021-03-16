#=============================================================================
#Combined plots for FVOI paper
#=============================================================================
load("fvoi_plot1.var")
library(gridExtra)
require(grid)
library(tidyverse)
#=============================================================================
# Paper plots
#=============================================================================
#=============================================================================
# Box 1
#=============================================================================
###Social info model
both_long_use = subset(both_long, time <= 60 )
both_long_use$species[both_long_use$species=="1"] = "species 1, no info"
both_long_use$species[both_long_use$species=="2"] = "species 2, no info"
both_long_use$species[both_long_use$species=="3"] = "species 1, social info"
both_long_use$species[both_long_use$species=="4"] = "species 2, social info"

#For text plotting
xpos = c(matrix(40,4,1))
ypos = c(both_long_use$N[both_long_use$time == 60])
ypos = ypos + (c(-1,1,-1,1))
suse = unique(both_long_use$species)

p0 =ggplot()+ geom_line( data = both_long_use, aes ( x = time, y = N, color = species)  )+ 
	geom_text( aes(x = xpos, y = ypos, label = suse, color = suse) ) +
	ylab("Population")+ xlab("Time")+   scale_colour_viridis_d()+
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)

p0
ggsave(file="fvoi_box1.pdf", p0)


#=============================================================================
# Box 2
#=============================================================================
###Draw species intrinsic fitness curve and plot them overtop the realized 
###env distribution: 
r1 = seq(0,1,length=ngens)
sp1 = exp(-0.5* ( (r1-env_fit$opt[1])/(env_fit$var[1]) )^2 )
sp2 = exp(-0.5* ( (r1-env_fit$opt[2])/(env_fit$var[1]) )^2 )
elt = data.frame( env_fit$env, r1, sp1, sp2) 
names(elt) = c("env", "r1","sp1","sp2")

p0=ggplot(data=elt)+geom_histogram(aes(x=env,y=..ncount..),color="black",fill="white")+
geom_line(aes(x=r1, y=sp1),color="#440154FF")+
geom_line(aes(x=r1, y=sp2),color="#35B779FF")+
ylab("Fitness value")+ xlab("Environment")+
theme_bw() + theme(
	text = element_text(size=10),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p0

###Environment, species fitness, germination response
ngens = dim(env_fit$fr)[1]
elt = data.frame( 1:ngens, env_fit$env, env_fit$fr, env_fit$gr,runif(ngens),runif(ngens)) 
names(elt) = c("Time", "env", "fr1","fr2","gr1","gr2","rgr1","rgr2")
#el_long = elt %>% gather(fr, repro, fr1:fr2) %>% gather(gr, germ, gr1:gr2)
#el_long = elt %>% gather(fr, repro, env:gr2) 
# el_long = elt %>% gather(fr, repro, fr1:rgr2) 
# el_long$fr[el_long$fr =="fr1" | el_long$fr =="gr1"|el_long$fr =="rgr1"] = "sp1"
# el_long$fr[el_long$fr =="fr2" | el_long$fr =="gr2"|el_long$fr =="rgr2"] = "sp2"
el_long = elt %>% gather(fr, repro, fr1:fr2) %>% gather(gr, germ, gr1:gr2)%>% 
gather(rgr, rgerm, rgr1:rgr2)
el_long$fr[el_long$fr =="fr1" ] = "sp1"
el_long$gr[el_long$gr =="gr1"] = "sp1"
el_long$rgr[el_long$rgr =="rgr1"] = "sp1"

el_long$fr[el_long$fr =="fr2" ] = "sp2"
el_long$gr[el_long$gr =="gr2"] = "sp2"
el_long$rgr[el_long$rgr =="rgr2"] = "sp2"

el2 = subset(el_long, Time<21)
c_use = c(color="#440154FF","#35B779FF" )

p1 = ggplot(data=el2) + geom_line(aes(x=Time, y=repro,color =fr ))+
geom_line(aes(x=Time, y=germ,color =gr ))+
geom_line(aes(x=Time, y=rgerm,color =rgr ),linetype = "dotted")+
scale_color_manual(values=c_use)+
geom_hline(yintercept=1 ,linetype = "dashed")+
 #scale_colour_viridis_d()+ 
 	ylab("Germination / Reproduction")+ 
	theme_bw() + theme(
	text = element_text(size=10),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p1

###Lottery model
ll_sub = subset(lott_long, time <= 60)
ll_sub$species[ll_sub$species=="1"] = "species 1, cue"
ll_sub$species[ll_sub$species=="2"] = "species 2, cue"
ll_sub$species[ll_sub$species=="3"] = "species 1, no cue"
ll_sub$species[ll_sub$species=="4"] = "species 2, no cue"


#For text plotting
xpos2 = c(matrix(10,4,1))
ypos2 = c(ll_sub$N[ll_sub$time == 1])
ypos2 = ypos2 + (c(0.2,-0.05,-0.1,0.5))
suse2 = unique(ll_sub$species)

p2 =ggplot()+ geom_line( data = ll_sub, aes ( x = time, y = N, color = species)  )+ 
	geom_text( aes(x = xpos2, y = ypos2, label = suse2, color = suse2) ) +
	ylab("Population")+ xlab("Time")+   scale_colour_viridis_d()+ 
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)

p2


# g=grid.arrange(p0, p1, p2, widths=c(unit(0.5, "npc"), unit(0.5, "npc") ),
# 					 heights=unit(0.5, "npc"), ncol = 2,
# 					 bottom = textGrob("Time",gp = gpar(fontsize = 14) ) )

# g= grid.arrange(arrangeGrob(p0,p1, ncol=1, nrow=2),
#          arrangeGrob(p2, ncol=1, nrow=1), heights=c(4,1), widths=c(1,2))

g= grid.arrange(arrangeGrob(p0,p1, ncol=1, nrow=2),
         arrangeGrob(p2, ncol=1, nrow=1), widths=c(unit(0.5, "npc"), 
         	unit(0.75, "npc") ), heights=unit(0.5, "npc"))


ggsave(file="fvoi_box2.pdf", g)

#=============================================================================
# Box 3
#=============================================================================


#=============================================================================
# Figure 2
#=============================================================================

#=============================================================================
#Base R plots: 
#=============================================================================

n_plot = 6000
col_use = c("black", "#440154FF")
col_use2 = c("#20A387FF","#FDE725FF")
plot(1:n_plot, cc_noi_out[1:n_plot,2], col = col_use[1], t="l",axes=F,cex.main=1.3,cex.lab=1.3, xlab = "Time", ylab="Population",ylim = c(0,40), xaxs="i",yaxs="i")
for (p in 2:(nPsp+1)){ 
	lines(1:n_plot,cc_noi_out[1:n_plot,p],col= col_use[p-1] )
	lines(1:n_plot,cc_i_out[1:n_plot,p],col= col_use2[p-1] )
}

axis(2, at=seq(0,40,10),cex.axis=1)
axis(1, at=seq(0,n_plot,2000),cex.axis=1)

#=============================================================================
# Figure X in ms? 
#=============================================================================
load("env_fit1.var")
rho_data = cbind(env_fit$mc2_all,env_fit$mr2_all,env_fit$mc3_all,env_fit$mr3_all)
rho_data = as.data.frame(rho_data)
colnames(rho_data) = c("species1_comp_info","species2_comp_info",
						"species1_comp_noinfo","species2_comp_noinfo",
						"species1_info","species2_info",
						"species1_noinfo","species2_noinfo")
#Convert to long: 
rho_long = rho_data %>% gather( species, rho )
p1 =ggplot()+ geom_line( data = rho_long, aes ( x = rho, y = N, color = species)  )+ 
	geom_text( aes(x = xpos2, y = ypos2, label = suse2, color = suse2) ) +
	ylab("")+ xlab("")+   scale_colour_viridis_d()+ 
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
