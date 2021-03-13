#=============================================================================
#Combined plots for FVOI paper
#=============================================================================
load("fvoi_plot1.var")
library(gridExtra)
require(grid)
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
	ylab("Population")+ xlab("")+   scale_colour_viridis_d()+
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)

p0

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


p1 =ggplot()+ geom_line( data = ll_sub, aes ( x = time, y = N, color = species)  )+ 
	geom_text( aes(x = xpos2, y = ypos2, label = suse2, color = suse2) ) +
	ylab("")+ xlab("")+   scale_colour_viridis_d()+ 
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)

p1

g=grid.arrange(p0, p1, widths=c(unit(0.5, "npc"), unit(0.5, "npc") ),
					 heights=unit(0.5, "npc"), ncol = 2,
					 bottom = textGrob("Time",gp = gpar(fontsize = 14) ) )

ggsave(file="fvoi_box1.pdf", g)

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
