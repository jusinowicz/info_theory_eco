#Fit models to see what predicts Biomass the best
#rDIT_eq = subset(rDIT, rMI >0.01 )
rDIT_eq = rDIT[c(-71,-68),]
rDIT_eq = subset(rDIT_eq, rMI >0.01 )

# rDIT_eq_con = rDIT_con[c(-71,-68),]
# rDIT_eq_con = subset(rDIT_eq_con, rMI >0.01 )


#Log-log 
l_rDIT_eq = cbind(log (rDIT_eq[,1:8]), rDIT_eq[,9])
lm_nspp = lm(l_rDIT_eq$nspp~l_rDIT_eq$Biomass)
lm_nspp_log = lm(l_rDIT_eq$nspp~l_rDIT_eq$Biomass)
lm_H=lm(l_rDIT_eq$shannon~l_rDIT_eq$Biomass)
lm_rS=lm(l_rDIT_eq$rS~l_rDIT_eq$Biomass)
#lm_rS=lm(l_rDIT_eq$Biomass~l_rDIT_eq$rS+I(l_rDIT_eq$rS^2))
lm_rCE=lm(l_rDIT_eq$rCE~l_rDIT_eq$Biomass)
lm_rMI=lm(l_rDIT_eq$rMI~l_rDIT_eq$Biomass)
lm_rCEMI=lm(l_rDIT_eq$rCE+l_rDIT_eq$rMI~l_rDIT_eq$Biomass)
lm_rSMI=lm(l_rDIT_eq$rS+l_rDIT_eq$rMI)

# gam_rS= gam (Biomass~s(rS),data=rDIT_eq )
# pr_gam_rS = predict(gam_rS, newdata =data.frame(rS = rDIT_eq$rS), type="response")
# #Predict data from models to fit to figure

pr_nspp_log = data.frame( Biomass =  exp(predict.lm ( lm_nspp_log) ),
	nspp = exp(l_rDIT_eq$nspp ) )
pr_H = data.frame( Biomass = exp(predict( lm_H) ) ,	shannon =exp(l_rDIT_eq$shannon ) )
pr_rS =data.frame( Biomass = exp(predict (lm_rS ) ), rS= exp(l_rDIT_eq$rS  ) )
pr_rCE = data.frame( Biomass = exp(predict(lm_rCE) ), rCE = exp(l_rDIT_eq$rCE ) )
pr_rMI = data.frame( Biomass = exp(predict(lm_rMI) ), rMI = exp(l_rDIT_eq$rMI ) )


summary(lm_nspp )
summary(lm_nspp_log )
summary(lm_H )
summary(lm_rS )
summary(lm_rCE )
summary(lm_rMI )
summary(lm_rCEMI )
summary(lm_rSMI )