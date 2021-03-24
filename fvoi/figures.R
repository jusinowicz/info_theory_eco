#=============================================================================
#Combined plots for FVOI paper
#=============================================================================
library(gridExtra)
require(grid)
library(tidyverse)
#=============================================================================
# Paper plots
#=============================================================================
#=============================================================================
# Box 1
#=============================================================================
load("fvoi_plot1.var")
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
load("ni_simple")
ngens = dim(Ni)[1]
ni = data.frame(1:ngens, Ni[,1], Ni_i[,1])
names(ni) = c("Time","ni1","ni_i1")
nil = ni%>%gather(ni, pop, ni1:ni_i1) #%>% gather(ni_i, pop_i, ni_i1:ni_i2)
c_use = c("#440154FF","#440154FF","#440154FF","#440154FF" )
ni2 = nil[nil$Time < 100,]

#For text plotting
xpos2 = c(matrix(c(55, 60), 2,1))
ypos2 = c(ni2$pop[ni2$Time == 40],ni2$pop[ni2$Time == 50])
ypos2= ypos2[c(1,4)]
suse2 = c("\u03C1(no information)", "\u03C1(information)")
#suse2 = c("A","B")

p1 = ggplot() + geom_line(data=ni2,aes(x=Time, y=pop,color =ni,linetype = ni ))+
	#geom_line(data=ni2,aes(x=Time, y=pop_i,color =ni_i ),linetype = "dotted")+
	scale_color_manual(values=c_use)+
	 #scale_colour_viridis_d()+ 
	ylab("Population")+ xlab("")+
	geom_text( aes(x = xpos2, y = ypos2, label = suse2,color=suse2)) +
 	scale_y_log10()+
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p1

load("dm_simp.var")
ngens = dim(Ni)[1]
ni = data.frame(1:ngens, Ni[,1], No[,1])
names(ni) = c("Time","ni1","no1")
nil = ni%>%gather(ni, pop, ni1:no1) #%>% gather(ni_i, pop_i, ni_i1:ni_i2)
c_use = c("#440154FF","#440154FF","#440154FF","#440154FF" )
ni2 = nil[nil$Time < 100,]

#For text plotting
xpos3 = c(matrix(c(55, 60), 2,1))
ypos3 = c(ni2$pop[ni2$Time == 40],ni2$pop[ni2$Time == 50])
ypos3= ypos3[c(1,4)]
suse3 = c( "\u03C1(information)","\u03C1(no information)")
#suse2 = c("A","B")

p2 = ggplot() + geom_line(data=ni2,aes(x=Time, y=pop,color =ni,linetype = ni ))+
	#geom_line(data=ni2,aes(x=Time, y=pop_i,color =ni_i ),linetype = "dotted")+
	scale_color_manual(values=c_use)+
	 #scale_colour_viridis_d()+ 
	ylab("")+ xlab("")+
	geom_text( aes(x = xpos3, y = ypos3, label = suse3,color=suse2)) +
 	scale_y_log10()+
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p2


g=grid.arrange(p1, p2, widths=c(unit(0.5, "npc"), unit(0.5, "npc") ),
					 heights=unit(0.5, "npc"), ncol = 2,
					 bottom = textGrob("Time",gp = gpar(fontsize = 14) ) )

ggsave(file="fvoi_box3.pdf",g)
#cairo_pdf(file="fvoi_box3.pdf")
#=============================================================================
# Figure 2
#=============================================================================
load("env_fit2.var")
ngens = dim(env_fit$mc2_all)[1]


rhos = NULL
for (s in 1:nspp){ 

	#Convert these to data frames
	mc21 = as.data.frame(env_fit$mc2_all[,,s])
	mr21 = as.data.frame(env_fit$mr2_all[,,s])
	mc31 = as.data.frame(env_fit$mc3_all[,,s])
	mr31 = as.data.frame(env_fit$mr3_all[,,s])

	#Gather all of the runs into long format
	mc21_long = mc21%>%gather(run1, val1,V1:V50) #Competition with info
	mr21_long = mr21%>%gather(run2, val2,V1:V50) #Competition without info
	mc31_long = mc31%>%gather(run3, val3,V1:V50) #No competition, info
	mr31_long = mr31%>%gather(run4, val4,V1:V50) #No competition, no info

	#Add a column to denote species
	gr_use = paste("s",s,sep="")
	mc21_long=cbind(mc21_long, gr = gr_use)
	# mr21_long=cbind(mr21_long, gr = gr_use)
	# mc31_long=cbind(mc31_long, gr = gr_use)
	# mr31_long=cbind(mr31_long, gr = gr_use)

	#Do the subtraction for the fvoi
	rhos_tmp = NULL
	rhos_tmp = cbind(mc21_long,mr21_long,mc31_long,mr31_long )
	rhos_tmp$val4 = rhos_tmp$val3 - rhos_tmp$val4
	rhos_tmp$val3 = rhos_tmp$val1 - rhos_tmp$val2

	niches = seq(1, 0, length = ngens) - 0.29
	rhos_tmp = cbind(Competition = niches, rhos_tmp)

	rhos = rbind(rhos,rhos_tmp)

}


r1 = rhos[rhos$Competition>niches[25],]
c_use = c(color="#440154FF","#35B779FF" )
#r1 = rhos

p0 = ggplot() +
	geom_point(data=r1, aes(x=Competition, y=val3, group = gr,color=gr ) )+
	geom_smooth(data=r1, method="lm" , formula = y ~ poly(x, 3), aes(x=Competition, y=val3, group = gr,color=gr ) )+
	geom_point(data=r1, aes(x=Competition, y=val4, group = gr,color=gr ) )+
	geom_smooth(data=r1, method="lm" , aes(x=Competition, y=val4, group = gr,color=gr ) )+
	scale_color_manual(values=c_use)+
	ylab("Fitness value of information")+ xlab("")+
	geom_hline(yintercept=0 ,linetype = "dashed")+
	#geom_text( aes(x = xpos2, y = ypos2, label = suse2,color=suse2)) +
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p0

p1 = ggplot() +
	geom_point(data=r1, aes(x=Competition, y=val1, group = gr,color=gr ) )+
	geom_smooth(data=r1, method="lm" , formula = y ~ poly(x, 3), aes(x=Competition, y=val1, group = gr,color=gr ) )+
	geom_point(data=r1, aes(x=Competition, y=val2, group = gr,color=gr ) )+
	geom_smooth(data=r1, method="lm" , formula = y ~ poly(x, 3),aes(x=Competition, y=val2, group = gr,color=gr ) )+
	scale_color_manual(values=c_use)+
	ylab("Fitness")+ xlab("")+
	geom_hline(yintercept=0 ,linetype = "dashed")+
	#geom_text( aes(x = xpos2, y = ypos2, label = suse2,color=suse2)) +
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p1

g=grid.arrange(p0, p1, widths=c(unit(0.5, "npc"), unit(0.5, "npc") ),
					 heights=unit(0.5, "npc"), ncol = 2,
					 bottom = textGrob("Niche overlap",gp = gpar(fontsize = 14) ) )

ggsave(file="fig2.pdf", g)

#=============================================================================
#Misc old plots: 
#=============================================================================
load("env_fit1.var")
rhos = data.frame(1:ngens,  env_fit$mc2_all, env_fit$mr2_all, env_fit$mc2_all-env_fit$mr2_all,env_fit$mc3_all-env_fit$mr3_all )
names(rhos) = c("Competition", "cc1","cc2", "cr1","cr2","crho1","crho2","rho1","rho2")
rhos_long = rhos %>% gather(cc, r_cc, cc1:cc2 )%>%
					gather(cr, r_cr, cr1:cr2 )%>%
					 gather(crho, r, crho1:crho2 )%>% 
					 gather(rho, rr, rho1:rho2 )

r1 = rhos_long[rhos_long$Competition<25,]
p0 = ggplot() +
geom_line(data=r1, aes(x=Competition, y=r, group = crho ) )+
geom_line(data=r1, aes(x=Competition, y=rr, group = rho ) )
p0


p1 = ggplot() + geom_line(data=r1, aes(x=Competition, y=r_cc, group = cc  ) )+
geom_line(data=r1, aes(x=Competition, y=r_cr, group = cr  ) )
p1

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


#Information gain as a function of model
plot(env_fit$mc2_all[,1]-env_fit$mr2_all[,1],col="red")                                                                  
points(env_fit$mc3_all[,1]-env_fit$mr3_all[,1])   

plot(env_fit$mc2_all[,2]-env_fit$mr2_all[,2],col="red")                                                                  
points(env_fit$mc3_all[,2]-env_fit$mr3_all[,2])   

#Gain in rho across competition.
plot(env_fit$mc2_all[,1],col="red", ylim=c(0,2) )                                                                  
points(env_fit$mr2_all[,1])       

plot(env_fit$mc2_all[,2],col="red", ylim=c(0,2) )                                                                  
points(env_fit$mr2_all[,2])       


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
