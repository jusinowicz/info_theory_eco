#=============================================================================
# R code for various information and dynamic information theoretic quantities. 
# These are meant to be used with multi-species population models representing
# modeled food webs.  
#
# The functions consist of two general types: 1) when the underlying 
# transition probabilities are known (as per the functions in 
# food_web_functions.R), and 2) when the data are time-series of population 
# counts. 
#
# The functions of type (2) discretize population counts and turn them into 
# blocks of symbols based on these counts which are treated as strings.  
#=============================================================================


#=============================================================================
# Time-series information theoretic functions
#=============================================================================
#=============================================================================
# Helper probability functions. These are functions to get various joint and 
# marginal probability distributions. These are used repeatedly in the various
# functions below. 
#=============================================================================
#=============================================================================
# get_k()
# Count blocks for block marginal probabilities. 
# series1			The full time series on which to do the calculations, as a
#					matrix or data frame. 
# k 				The size of the block history to be used. 
# width_A			The largest width (i.e. maximum number of digits) of any
#					word.  
#=============================================================================

get_k = function (series1,k,width_A){

	totk = dim(as.matrix(series1))[1] - k+1
	blocks_k = matrix(0.0, totk, 1)

	for (n in 1:totk) { 
		#Use width_A to format the "names" of blocks by padding zeros as needed
		blocks_k[n] = paste(c (formatC(series1[n:((n+k-1))],width=width_A, 
			flag = "0" )),collapse="")
	}

	return(blocks_k)
}

#=============================================================================
# get_kpf()
# Count the joint historical and future blocks. 
# series1			The full time series on which to do the calculations, as a
#					matrix or data frame. 
# k 				The size of the block history to be used.
# s 				The size of the block future to be used (defaults to k)
# width_A			The largest width (i.e. maximum number of digits) of any
#					word.  
#=============================================================================

get_kpf = function (series1, k=2, s=k, width_A) {
	k2 = k+s
	totk = dim(as.matrix(series1))[1] - k2+1
	blocks_kpf = matrix(0.0, totk, 1) #joint of X(k)+1, X(k), Y(l)

	for (n in 1:totk) { 
		#In order, this is X(s),X(k)
		blocks_kpf[n] = paste(c (formatC(series1[ (n+k):(n+k2-1)],width=width_A, 
			flag = "0" ), "i", 
		formatC(series1[(n:(n+k-1))],width=width_A, flag = "0" ) ),collapse="")

	}

	# for (n in 1:totk) { 
	# 	#In order, this is X(s),X(k)
	# 	blocks_kpf[n] = paste(c (formatC (series1[(n:(n+k-1))],width=width_A, flag = "0" ) ,
	# 	 "i", formatC(series1[ (n+k):(n+k2-1)],width=width_A, flag = "0" )
	# 	 ),collapse="")

	# }

	return(blocks_kpf)

}


#=============================================================================
# get_kY()
# Count neighborhood-history blocks, mainly for the Transfer Entropy. This is 
# is a joint block over (X(k), Y(l) ) 
#	X(k) is the historical block
#	Y(l) is the value of all neighbors at X(k)
#	
# series1			The full time series on which to do the calculations, as a
#					matrix or data frame. 
# k 				The size of the block history to be used.
# focal				Which column in series1 is the time series of interest? 
#					The default is to assume it is the first column.
# width_A			The largest width (i.e. maximum number of digits) of any
#					word.  
#=============================================================================

get_kY = function (series1, k=2, focal=1, width_A) {

	totk = dim(series1)[1] - k+1
	blocks_kY = matrix(0.0, totk, 1) #joint of X(k)+1, X(k), Y(l)
	Yl=(1:dim(series1)[2])[-focal]

	for (n in 1:totk) { 
		#In order, this is X(k),Y(l)
		blocks_kY[n] =  paste(c (formatC(series1[(n:(n+k-1)),focal],width=width_A, 
			flag = "0" ),"i", 
		formatC(series1[(n+k-1),Yl],width=width_A, flag = "0" )  ),collapse="")

	}

	return(blocks_kY)

}

#=============================================================================
# get_1kY()
# Count blocks for the Transfer Entropy. This is is a joint block over 
# the triplet (X(k)+1, X(k), Y(l) ) 
#	X(k) is the historical block
# 	X(k)+1 is the next value after the block
#	Y(l) is the value of all neighbors at X(k)
#	
# series1			The full time series on which to do the calculations, as a
#					matrix or data frame. 
# k 				The size of the block history to be used.
# focal				Which column in series1 is the time series of interest? 
#					The default is to assume it is the first column.
# width_A			The largest width (i.e. maximum number of digits) of any
#					word.  
#=============================================================================

get_1kY = function (series1, k=2, focal=1, width_A) {
	k2 = k+1
	totk = dim(series1)[1] - k2+1
	blocks_1kY = matrix(0.0, totk, 1) #joint of X(k)+1, X(k), Y(l)
	Yl=(1:dim(series1)[2])[-focal]

	for (n in 1:totk) { 
		#In order, this is X(k)+1, X(k),Y(l)
		blocks_1kY[n] =  paste(c (formatC(series1[n+k2-1,focal],width=width_A, 
			flag = "0" ), "i", 
		formatC(series1[(n:(n+k-1)),focal],width=width_A, flag = "0" ),"i",
		formatC(series1[(n+k-1),Yl],width=width_A, flag = "0" )  ), 
		collapse="")

	}

	return(blocks_1kY)

}

#=============================================================================
# Dynamic information quantities! 
#=============================================================================
#=============================================================================
# get_MI()
# Mutual Information
#
#=============================================================================


#=============================================================================
# get_EE()
# Excess entropy
# series1			The full time series on which to do the calculations, as a
#					matrix or data frame. 
# k 				The size of the block history to be used. 
# focal				Which column in series1 is the time series of interest? 
#					The default is to assume it is the first column. 
# with_blocks=F		When set to TRUE, the function will also append all of 
#					the lists of blocks.
# blocks_only=F 	When set to TRUE, the function will ONLY return the 
#					lists of blocks. 
#=============================================================================
get_EE = function ( series1, k=2, s=k, focal = 1, with_blocks=FALSE, 
	blocks_only=FALSE){
	
	EE1 = NULL #Attach both the average and the local version to this

	k2=k+s
	ngensa = dim(series1)[1] - k+1
	ngensb = dim(series1)[1] - k2+1
	blocks_k = matrix(0.0, ngensa, 1) # marg of X(k)
	blocks_kpf = matrix(0.0, ngensb, 1) #joint of X(k)+s, X(k) 

	#Need some additional pre-formatting to make this work: 
	#Find number of digits in largest integer: 
	width_A = nchar(trunc(max(d1_use)))

	######################
	# p( X(k) )
	######################
	blocks_k = get_k(series1[,focal],k,width_A )
	marg_k =  ( prop.table(table( as.data.frame(blocks_k) )))

	######################
	# p( X(k)+s,X(k) )
	######################
	blocks_kpf = get_kpf(series1[,focal],k,s,width_A)
	joint_kpf = ( prop.table(table( as.data.frame(blocks_kpf) )))

	#Calculate Excess Entropy
	nblocks = length(joint_kpf)
	EE1_mean =matrix(0.0, nblocks, 1)
	for (n in 1:nblocks){ 

		#Split the whole block into its subblock, X(k)+1, X(k), and Y(l)
		#The right/future block
		block1 = unlist(strsplit(names(joint_kpf)[n],"i") )[1] 
		#The left/past block
		block2 = unlist(strsplit(names(joint_kpf)[n],"i") )[2]

		EE1_mean[n] = joint_kpf[n]*log2(joint_kpf[n]/( 
			marg_k[block2]*marg_k[block1]  )) 
	}

	#Calculate the local Excess Entropy
	nblocks = length(blocks_kpf)
	EE1_local =matrix(0.0, nblocks, 1)
	for (n in 1:nblocks){ 

		#Split the whole block into its subblock, X(k)+s, X(k)
		#The right/future block
		block1 = unlist(strsplit(blocks_kpf[n],"i") )[1] 
		#The left/past block
		block2 = unlist(strsplit(blocks_kpf[n],"i") )[2]

		EE1_local[n] = log2(joint_kpf[blocks_kpf[n]]/( 
			marg_k[block2]*marg_k[block1]  )) 
	}

	if (with_blocks == TRUE ){
		EE1$mean = sum(EE1_mean)
		EE1$local = EE1_local
		EE1$blocks_k = blocks_k
		EE1$blocks_kpf = blocks_kpf


	} else if (blocks_only == TRUE ){
		EE1$blocks_k = blocks_k
		EE1$blocks_kpf = blocks_kpf
	} else {  
		EE1$mean = sum(EE1_mean)
		EE1$local = EE1_local
	}
	
	return(EE1)

}

#=============================================================================
# get_ais()
# Active information storage
# series1			The full time series on which to do the calculations, as a
#					matrix or data frame. 
# k 				The size of the block history to be used. 
# focal				Which column in series1 is the time series of interest? 
#					The default is to assume it is the first column. 
# with_blocks=F		When set to TRUE, the function will also append all of 
#					the lists of blocks.
# blocks_only=F 	When set to TRUE, the function will ONLY return the 
#					lists of blocks. 
#=============================================================================
get_ais = function ( series1, k=2, s=1, focal = 1, with_blocks=FALSE, 
	blocks_only=FALSE){
	
	AIS1 = NULL #Attach both the average and the local version to this

	k2=k+s 
	ngensa = dim(series1)[1] - k+1
	ngensb = dim(series1)[1] - k2+1
	blocks_1 = matrix(0.0, ngensb, 1) #marg of X(k)+1
	blocks_k = matrix(0.0, ngensa, 1) # marg of X(k)
	blocks_1k = matrix(0.0, ngensb, 1) #joint of X(k)+1, X(k) 

	#Need some additional pre-formatting to make this work: 
	#Find number of digits in largest integer: 
	width_A = nchar(trunc(max(d1_use)))

	######################
	# p( X(k)+1 )
	######################
	blocks_1 = get_k(series1[(k2):dim(series1)[1],focal],k=1,width_A)
	marg_1 =  ( prop.table(table( as.data.frame(blocks_1) )))

	######################
	# p( X(k) )
	######################
	blocks_k = get_k(series1[,focal],k,width_A )
	marg_k =  ( prop.table(table( as.data.frame(blocks_k) )))

	######################
	# p( X(k)+1,X(k) )
	######################
	blocks_1k = get_kpf(series1[,focal],k,s,width_A)
	joint_1k = ( prop.table(table( as.data.frame(blocks_1k) )))

	#Calculate active information storage
	nblocks = length(joint_1k)
	AIS1_mean =matrix(0.0, nblocks, 1)
	for (n in 1:nblocks){ 

		#Split the whole block into its subblock, X(k)+1, X(k), and Y(l)
		#The right/future block
		block1 = unlist(strsplit(names(joint_1k)[n],"i") )[1] 
		#The left/past block
		block2 = unlist(strsplit(names(joint_1k)[n],"i") )[2]

		AIS1_mean[n] = joint_1k[n]*log2( (joint_1k[n])/( 
			marg_k[ block2 ]*marg_1[ block1 ] ) ) 
		#AIS1[n] = log2( (joint_1k[n])/( marg_k[ block2 ]*marg_1[ block1 ] ) ) 

	}

	#Calculate the local AIS
	nblocks = length(blocks_1k)
	AIS1_local =matrix(0.0, nblocks, 1)
	for (n in 1:nblocks){ 

		#Split the whole block into its subblock, X(k)+1, X(k), and Y(l)
		#The right/future block
		block1 = unlist(strsplit(blocks_1k[n],"i") )[1] 
		#The left/past block
		block2 = unlist(strsplit(blocks_1k[n],"i") )[2]

		AIS1_local[n] = log2( (joint_1k[blocks_1k[n]])/( 
			marg_k[ block2 ]*marg_1[ block1 ] ) ) 

	}

	if (with_blocks == TRUE ){
		AIS1$mean = sum(AIS1_mean)
		AIS1$local = AIS1_local
		AIS1$blocks_1=blocks_1
		AIS1$blocks_k = blocks_k
		AIS1$blocks_1k = blocks_1k

	} else if (blocks_only == TRUE ){
		AIS1$blocks_1=blocks_1
		AIS1$blocks_k = blocks_k
		AIS1$blocks_1k = blocks_1k
	} else {  
		AIS1$mean = sum(AIS1_mean)
		AIS1$local = AIS1_local
	}
	
	return(AIS1)


}


#=============================================================================
# get_TE()
# The Transfer entropy. Which version of the TE gets calculated depends on 
# what gets passed to the function. E.g. when all other linked nodes, agents, 
# etc. are given in series1 then this is the complete transfer entropy. 
# Otherwise it will be some version of the apparent TE.   
# series1			The full time series on which to do the calculations, as a
#					matrix or data frame. 
# k 				The size of the block history to be used. 
# focal				Which column in series1 is the time series of interest? 
#					The default is to assume it is the first column. 
# with_blocks=F		When set to TRUE, the function will also append all of 
#					the lists of blocks.
# blocks_only=F 	When set to TRUE, the function will ONLY return the 
#					lists of blocks. 
#=============================================================================

get_TE = function ( series1, k=2, s=1, focal = 1, with_blocks=FALSE, 
	blocks_only=FALSE){
	
	TE1 = NULL #Attach both the average and the local version to this

	#Make l everything but the source series
	Yl=(1:dim(series1)[2])[-focal]
	k2=k+1 
	ngensa = dim(series1)[1] - k+1
	ngensb = dim(series1)[1] - k2+1
	blocks_k = matrix(0.0, ngensa, 1) #Block marginals
	blocks_1kY = matrix(0.0, ngensb, 1) #joint of X(k)+1, X(k), Y(l)
	blocks_kY = matrix(0.0, ngensa, 1) #joint of X(k), Y(l)
	blocks_1k = matrix(0.0, ngensb, 1) #joint of X(k)+1, X(k)

	#Need some additional pre-formatting to make this work: 
	#Find number of digits in largest integer: 
	width_A = nchar(trunc(max(d1_use)))

	#There are four total probabilities that are calculated here. 
	#These include three different joint probability distributions 
	#and one marginal distribution
	######################
	# p( X(k) )
	######################
	blocks_k = get_k(series1[,focal],k,width_A )
	marg_k =  ( prop.table(table( as.data.frame(blocks_k) )))

	######################
	# p( X(k)+1,X(k),Y(l) )
	######################
	blocks_1kY = get_1kY(series1,k,focal,width_A ) 
	joint_1kY =   ( prop.table(table( as.data.frame(blocks_1kY) )))

	######################
	# p( X(k),Y(l) )
	######################
	blocks_kY = get_kY(series1,k,focal,width_A ) 
	joint_kY =  ( prop.table(table( as.data.frame(blocks_kY) )))

	######################
	# p( X(k)+1,X(k) )
	######################
	blocks_1k = get_kpf(series1[,focal],k,s,width_A)
	joint_1k = ( prop.table(table( as.data.frame(blocks_1k) )))

	#Now calculate the transfer entropy
	nblocks = length(joint_1kY)
	TE1_mean =matrix(0.0, nblocks, 1)
	for (n in 1:nblocks){ 

		#Split the whole block into its subblock, X(k)+1, X(k), and Y(l)
		#The right/future block
		block1 = unlist(strsplit(names(joint_1kY)[n],"i") )[1] 
		#The left/past block
		block2 = unlist(strsplit(names(joint_1kY)[n],"i") )[2]
		#The neighborhood block
		block3 = unlist(strsplit(names(joint_1kY)[n],"i") )[3]

		TE1_mean[n] = joint_1kY[n]*log2( (joint_1kY[n]*marg_k[block2])/( 
			joint_kY[ paste( c(block2,"i",block3) ,collapse="") ]*
			joint_1k[paste( c(block1,"i",block2) ,collapse="") ] ) ) 

	}	

	#Calculate the local transfer entropy
	nblocks = length(blocks_1kY)
	TE1_local =matrix(0.0, nblocks, 1)
	for (n in 1:nblocks){ 

		#Split the whole block into its subblock, X(k)+1, X(k), and Y(l)
		#The right/future block
		block1 = unlist(strsplit(blocks_1kY[n],"i") )[1] 
		#The left/past block
		block2 = unlist(strsplit(blocks_1kY[n],"i") )[2]
		#The neighborhood block
		block3 = unlist(strsplit(blocks_1kY[n],"i") )[3]

		TE1_local[n] = log2( (joint_1kY[blocks_1kY[n]]*marg_k[block2])/( 
			joint_kY[ paste( c(block2,"i",block3) ,collapse="") ]*
			joint_1k[paste( c(block1,"i",block2) ,collapse="") ] ) ) 

	}	

	if (with_blocks == TRUE ){
		TE1$mean = sum(TE1_mean)
		TE1$local = TE1_local
		TE1$blocks_k=blocks_k
		TE1$blocks_kY = blocks_kY
		TE1$blocks_1k = blocks_1k
		TE1$blocks_1kY = blocks_1kY

	} else if (blocks_only == TRUE ){
		TE1$blocks_k=blocks_k
		TE1$blocks_kY = blocks_kY
		TE1$blocks_1k = blocks_1k
		TE1$blocks_1kY = blocks_1kY

	} else {  
		TE1$mean = sum(TE1_mean)
		TE1$local = TE1_local
	}
	
 	return( TE1 )
	

}

#=============================================================================
# get_SI()
# The Separable Information 
# series1			The full time series on which to do the calculations, as a
#					matrix or data frame. 
# k 				The size of the block history to be used. 
# focal				Which column in series1 is the time series of interest? 
#					The default is to assume it is the first column. 
# with_blocks=F		When set to TRUE, the function will also append all of 
#					the probability disributions of blocks.
# blocks_only=F 	When set to TRUE, the function will ONLY return the 
#					probability distributions of blocks. 
#=============================================================================

get_SI= function ( series1, k=2, s=1, focal = 1, with_blocks=FALSE, 
	blocks_only=FALSE){
	
	SI1 = NULL #Attach both the average and the local version to this
	
	#Make l everything but the source series
	Yl=(1:dim(series1)[2])[-focal]
	
	#First calculate the active information storage:
	AIS_temp = get_ais ( series1, k=k, s=s, focal = focal, with_blocks=with_blocks,
		blocks_only=blocks_only )
	
	#Sum over the apparent transfer entropy between series1 and each other node
	TE_temp = NULL
	TE_temp$mean = NULL
	TE_temp$local = NULL
	blocks_k = list(matrix(0.0, Yl, 1)) #Block marginals
	blocks_1kY = list(matrix(0.0, Yl, 1)) #joint of X(k)+1, X(k), Y(l)
	blocks_kY = list(matrix(0.0, Yl, 1)) #joint of X(k), Y(l)
	blocks_1k = list(matrix(0.0, Yl, 1)) #joint of X(k)+1, X(k)
 

	for (i in 1:length(Yl)){
		TE_temp2 = get_TE( series1[,c(focal,Yl[i])],k=k,s=s,focal=1,
			with_blocks=with_blocks, blocks_only=blocks_only )
		
		if (with_blocks == TRUE ){

			if(is.null(TE_temp$mean)){
				TE_temp$mean = TE_temp2$mean
				TE_temp$local = TE_temp2$local
			} else {
				TE_temp$mean = TE_temp$mean + TE_temp2$mean
				TE_temp$local = TE_temp$local + TE_temp2$local
			}

			blocks_k[i]=list(blocks_k)
			blocks_kY[i] = list(blocks_kY)
			blocks_1k[i] = list(blocks_1k)
			blocks_1kY[i] = list(blocks_1kY)

		} else if (blocks_only == TRUE ){
			blocks_k[i]=list(blocks_k)
			blocks_kY[i] = list(blocks_kY)
			blocks_1k[i] = list(blocks_1k)
			blocks_1kY[i] = list(blocks_1kY)

		} else {  

			if(is.null(TE_temp$mean)){
				TE_temp$mean = TE_temp2$mean
				TE_temp$local = TE_temp2$local
			} else {
				TE_temp$mean = TE_temp$mean + TE_temp2$mean
				TE_temp$local = TE_temp$local + TE_temp2$local
			}	
		}


		
		

	} 
	
	if (with_blocks == TRUE ){
		SI1$mean = AIS_temp$mean + TE_temp$mean
		SI1$local = AIS_temp$local + TE_temp$local
		SI1$blocks_k=blocks_k
		SI1$blocks_kY = blocks_kY
		SI1$blocks_1k = blocks_1k
		SI1$blocks_1kY = blocks_1kY

	} else if (blocks_only == TRUE ){
		SI1$blocks_k=blocks_k
		SI1$blocks_kY = blocks_kY
		SI1$blocks_1k = blocks_1k
		SI1$blocks_1kY = blocks_1kY

	} else {  
		SI1$mean = AIS_temp$mean + TE_temp$mean
		SI1$local = AIS_temp$local + TE_temp$local
	}

 	return( SI1 )
	 

}


#=============================================================================
#get_info_dynamics
# This is mostly a wrapper to apply and sort the dynamic info metrics to a time
# series of an entire foodweb.
#
# pop_ts				Population matrix with time as rows, each species as column
# k 				The size of the block history to be used. 
# with_blocks=F		When set to TRUE, the function will also append all of 
#					the probability disributions of blocks.
# blocks_only=F 	When set to TRUE, the function will ONLY return the 
#					probability distributions of blocks. 
#=============================================================================

get_info_dynamics = function (pop_ts, k=2,with_blocks=FALSE, blocks_only=FALSE ) { 

	d1_web = NULL

	ngens = dim(pop_ts)[1] #Time
	nspp = dim(pop_ts)[2] #Number of species

	if (with_blocks == TRUE ){
	
		#Use these to store the local and global average values
		di_ee_means = matrix(0,nspp,1)
		di_ai_means = matrix(0,nspp,1)
		di_te_means = matrix(0,nspp,1)
		di_si_means = matrix(0,nspp,1)

		
		di_ee_local = matrix(0,(ngens-2*k+1),nspp)
		di_ai_local = matrix(0,(ngens-k),nspp)
		di_te_local = matrix(0,(ngens-k),nspp)
		di_si_local = matrix(0,(ngens-k),nspp)

		di_ee_blocks = list(matrix(0,nspp,1))
		di_ai_blocks = list(matrix(0,nspp,1))
		di_te_blocks= list(matrix(0,nspp,1))
		di_si_blocks = list(matrix(0,nspp,1))

		for ( n in 1:nspp) {
			print(n)
			di_ee_temp = get_EE (pop_ts,k=k,focal=n,with_blocks=with_blocks, blocks_only=blocks_only)
			di_ai_temp = get_ais (pop_ts,k=k,focal=n,with_blocks=with_blocks, blocks_only=blocks_only)
			di_te_temp = get_TE (pop_ts,k=k,focal=n,with_blocks=with_blocks, blocks_only=blocks_only)
			di_si_temp = get_SI (pop_ts,k=k,focal=n,with_blocks=with_blocks, blocks_only=blocks_only)

			di_ee_means[n] = di_ee_temp$mean
			di_ai_means[n] = di_ai_temp$mean
			di_te_means[n] = di_te_temp$mean
			di_si_means[n] = di_si_temp$mean
			
			di_ee_local[,n] = di_ee_temp$local
			di_ai_local[,n] = di_ai_temp$local
			di_te_local[,n] = di_te_temp$local
			di_si_local[,n] = di_si_temp$local

			di_ee_blocks[n] = list(di_ee_temp[grepl("block",names(di_ee_temp))])
			di_ai_blocks[n] = list(di_ai_temp[grepl("block",names(di_ai_temp))])
			di_te_blocks[n] = list(di_te_temp[grepl("block",names(di_te_temp))])
			di_si_blocks[n] = list(di_si_temp[grepl("block",names(di_si_temp))])


		}

		#Mean quantities
		d1_web$ee_means = di_ee_means
		d1_web$ai_means = di_ai_means
		d1_web$te_means = di_te_means
		d1_web$si_means = di_si_means

		#Local quantities
		d1_web$ee_local = di_ee_local
		d1_web$ai_local = di_ai_local
		d1_web$te_local = di_te_local
		d1_web$si_local = di_si_local


		d1_web$ee_blocks = di_ee_blocks 
		d1_web$ai_blocks = di_ai_blocks 
		d1_web$te_blocks = di_te_blocks 
		d1_web$si_blocks = di_si_blocks 


	} else if (blocks_only == TRUE ){

		di_ee_blocks = list(matrix(0,nspp,1))
		di_ai_blocks = list(matrix(0,nspp,1))
		di_te_blocks= list(matrix(0,nspp,1))
		di_si_blocks = list(matrix(0,nspp,1))

		for ( n in 1:nspp) {
			print(n)
			di_ee_temp = get_EE (pop_ts,k=k,focal=n,with_blocks=with_blocks, blocks_only=blocks_only)
			di_ai_temp = get_ais (pop_ts,k=k,focal=n,with_blocks=with_blocks, blocks_only=blocks_only)
			di_te_temp = get_TE (pop_ts,k=k,focal=n,with_blocks=with_blocks, blocks_only=blocks_only)
			di_si_temp = get_SI (pop_ts,k=k,focal=n,with_blocks=with_blocks, blocks_only=blocks_only)

			di_ee_blocks[n] = list(di_ee_temp[grepl("block",names(di_ee_temp))])
			di_ai_blocks[n] = list(di_ai_temp[grepl("block",names(di_ai_temp))])
			di_te_blocks[n] = list(di_te_temp[grepl("block",names(di_te_temp))])
			di_si_blocks[n] = list(di_si_temp[grepl("block",names(di_si_temp))])


		}

		d1_web$ee_blocks = di_ee_blocks 
		d1_web$ai_blocks = di_ai_blocks 
		d1_web$te_blocks = di_te_blocks 
		d1_web$si_blocks = di_si_blocks 


	} else { 
	#Use these to store the local and global average values
		di_ee_means = matrix(0,nspp,1)
		di_ai_means = matrix(0,nspp,1)
		di_te_means = matrix(0,nspp,1)
		di_si_means = matrix(0,nspp,1)

		
		di_ee_local = matrix(0,(ngens-2*k+1),nspp)
		di_ai_local = matrix(0,(ngens-k),nspp)
		di_te_local = matrix(0,(ngens-k),nspp)
		di_si_local = matrix(0,(ngens-k),nspp)


		for ( n in 1:nspp) {
			print(n)
			di_ee_temp = get_EE (pop_ts,k=k,focal=n,with_blocks=with_blocks, blocks_only=blocks_only)
			di_ai_temp = get_ais (pop_ts,k=k,focal=n,with_blocks=with_blocks, blocks_only=blocks_only)
			di_te_temp = get_TE (pop_ts,k=k,focal=n,with_blocks=with_blocks, blocks_only=blocks_only)
			di_si_temp = get_SI (pop_ts,k=k,focal=n,with_blocks=with_blocks, blocks_only=blocks_only)

			di_ee_means[n] = di_ee_temp$mean
			di_ai_means[n] = di_ai_temp$mean
			di_te_means[n] = di_te_temp$mean
			di_si_means[n] = di_si_temp$mean
			
			di_ee_local[,n] = di_ee_temp$local
			di_ai_local[,n] = di_ai_temp$local
			di_te_local[,n] = di_te_temp$local
			di_si_local[,n] = di_si_temp$local


		}

		#Mean quantities
		d1_web$ee_means = di_ee_means
		d1_web$ai_means = di_ai_means
		d1_web$te_means = di_te_means
		d1_web$si_means = di_si_means

		#Local quantities
		d1_web$ee_local = di_ee_local
		d1_web$ai_local = di_ai_local
		d1_web$te_local = di_te_local
		d1_web$si_local = di_si_local

	}
	

	return(d1_web)
}

#=============================================================================
#get_te_web
# This is a wrapper to cycle through and get the Transfer Entropy on a pairwise
# basis for all species pairs in the food web. 
#
# pop_ts			Population matrix with time as rows, each species as column
# k 				The size of the block history to be used. 
#=============================================================================

get_te_web = function (pop_ts, k=2, s=1 ) { 

	ngens = dim(pop_ts)[1] #Time
	nspp = dim(pop_ts)[2] #Number of species
	
	te_web = matrix(0,nspp,nspp)


	for ( n in 1:nspp) {
		for ( n2 in 1:nspp) {
			print( paste( n," ",n2) )
			te_temp = get_TE( pop_ts[,c(n,n2)],k=k,s=s,focal=1)
			te_web[n,n2] = te_temp$mean
		}
	}
	return(te_web)
}

#=============================================================================
#get_si_web
# This is a wrapper to cycle through and get the Separable Information on a 
# pairwise basis for all species pairs in the food web. 
#
# pop_ts			Population matrix with time as rows, each species as column
# k 				The size of the block history to be used. 
#=============================================================================

get_si_web = function (pop_ts, k=2, s=1 ) { 

	ngens = dim(pop_ts)[1] #Time
	nspp = dim(pop_ts)[2] #Number of species
	
	si_web = matrix(0,nspp,nspp)


	for ( n in 1:nspp) {
		for ( n2 in 1:nspp) {
			print( paste( n," ",n2) )
			si_temp = get_SI( pop_ts[,c(n,n2)],k=k,s=s,focal=1)
			si_web[n,n2] = si_temp$mean
		}
	}
	return(si_web)
}