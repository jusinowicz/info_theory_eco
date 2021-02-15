#=============================================================================
#Combined plots for FVOI paper
#=============================================================================
#=============================================================================
# Paper plots
#=============================================================================
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

