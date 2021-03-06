#=============================================================================
# A series of functions to retrieve data from food web databases, including
# EcoBase and Mangal. 
# 
# Credit goes to Frederico Mestre:
# https://www.r-bloggers.com/function-to-download-biotic-interaction-datasets/
#=D============================================================================
library(RCurl)
library(XML)
library(plyr)
library(rmangal)
#=============================================================================
#get_eb 
#
# This function gets all of the complete, publicly available food webs from 
# the EcoBase database. It returns an object that includes each food web as 
# a member in a list. Each food web itself is a list of two objects by default:
# a list of each species' biomass and the interaction (consumption) matrix 
# of proportion consumed. 
#=============================================================================

#EcoBase
get_eb = function (biomass = TRUE, pb = TRUE, qb = TRUE, ee = TRUE, gs = TRUE, export= TRUE, 
    study_id = FALSE, ecosyst = FALSE, ref = FALSE ) {

  fwlist <- list()
    
    message("####################### ECOBASE DATABASE #######################\n\n")
    
    message("Fetching info from the EcoBase website!")
    suppressWarnings({
      
      #To obtain the list of available models
      
      suppressMessages({
        h=basicTextGatherer()
        curlPerform(url = 'http://sirs.agrocampus-ouest.fr/EcoBase/php/webser/soap-client_3.php',writefunction=h$update)
        data1 <- xmlTreeParse(h$value(),useInternalNodes=TRUE)
        liste_mod <- ldply(xmlToList(data1),data.frame)#liste_mod contains a list and decription
      })
      
      #Select only those allowing dissemination
      l2 <- subset(liste_mod, model.dissemination_allow =="true")#only those of which dissemination is allowed
      message("Sellected only those to which model dissemination is allowed!")
      
      #Select only those with whole food webs
      l3 <- subset(l2, model.whole_food_web =="true")#only those with the full food web
      message("Sellected only those to which the whole food web is available!")
      
      #Get model names
      model.name <- as.character(l3$model.model_name)
      input_list <- list()
      id <- as.numeric(as.character(l3$model.model_number))
      
      #Loop to get input list
      for(i in 1:nrow(l3)){
        
        message(paste0("Fetching information on food web ",i, " of ", nrow(l3)))
        
        suppressMessages({
          h=basicTextGatherer()
          mymodel <- id[i]
          curlPerform(url = paste('http://sirs.agrocampus-ouest.fr/EcoBase/php/webser/soap-client.php?no_model=',mymodel,sep=''),writefunction=h$update,verbose=TRUE)
          data2 <- xmlTreeParse(h$value(),useInternalNodes=TRUE)
          input1 <- xpathSApply(data2,'//group',function(x) xmlToList(x))
        })
        
        #need do name the columns
        names_input <- as.character(input1[1,])
        input1 <- as.data.frame(input1)
        colnames(input1) <- names_input
        input1 <- input1[-1,]
        input_list[[i]] <- input1
        print(i)
        print(input1[[1]]$habitat_area)
      }#end of loop to get input list
      
      mnames <- names(input_list)
      
      for (i in 1:length(input_list)){
        print(i)
        m2 <- input_list[[i]] #get the model
        
        nnodes <- length(m2)
        node_names <- names(m2)
        
        if (biomass == TRUE)
           { 
           nodes_biomass <- as.data.frame(matrix(ncol=3, nrow=nnodes))
           names(nodes_biomass) <- c("id", "name", "biomass")
           }

        if (pb == TRUE)
           { 
           nodes_pb <- as.data.frame(matrix(ncol=3, nrow=nnodes))
           names(nodes_pb) <- c("id", "name", "pb")
           }

        if (qb == TRUE)
           { 
           nodes_qb <- as.data.frame(matrix(ncol=3, nrow=nnodes))
           names(nodes_qb) <- c("id", "name", "qb")
           }
                 
        if (ee == TRUE)
           { 
           nodes_ee <- as.data.frame(matrix(ncol=3, nrow=nnodes))
           names(nodes_ee) <- c("id", "name", "ee")
           }

        if (gs == TRUE)
           { 
           nodes_gs <- as.data.frame(matrix(ncol=3, nrow=nnodes))
           names(nodes_gs) <- c("id", "name", "gs")
           }
                         
        if (export == TRUE)
           { 
           nodes_ex <- as.data.frame(matrix(ncol=3, nrow=nnodes))
           names(nodes_ex) <- c("id", "name", "ex")
           }
        
        
        
        int_matrix <- as.data.frame(matrix(ncol=nnodes, nrow=nnodes))
        
        for(j in 1:length(m2)){
          
          node1 <- m2[[j]]
          node_id <- as.numeric(node1$group_seq)  
          node1_biomass <- as.numeric(node1$biomass)/as.numeric(node1$habitat_area) 
          node1_pb <- as.numeric(node1$pb)  
          node1_qb <- as.numeric(node1$qb)  
          node1_ee <- as.numeric(node1$ee)  
          node1_gs <- as.numeric(node1$gs)  
          node1_ex <- as.numeric(node1$export)  


          node_name <- node_names[j]
          
          #biomass
          if (biomass == TRUE)
          {  
          nodes_biomass[node_id, 1] <- node_id  
          nodes_biomass[node_id, 2] <- node_name  
          nodes_biomass[node_id, 3] <- node1_biomass  
          }

          if (pb == TRUE)
          {  
          nodes_pb[node_id, 1] <- node_id  
          nodes_pb[node_id, 2] <- node_name  
          nodes_pb[node_id, 3] <- node1_pb
          }
          
          if (qb == TRUE)
          {  
          nodes_qb[node_id, 1] <- node_id  
          nodes_qb[node_id, 2] <- node_name  
          nodes_qb[node_id, 3] <- node1_qb 
          }
          
          if (ee == TRUE)
          {  
          nodes_ee[node_id, 1] <- node_id  
          nodes_ee[node_id, 2] <- node_name  
          nodes_ee[node_id, 3] <- node1_ee 
          }

         if (gs == TRUE)
          {  
          nodes_gs[node_id, 1] <- node_id  
          nodes_gs[node_id, 2] <- node_name  
          nodes_gs[node_id, 3] <- node1_gs 
          }
          
         if (export == TRUE)
          {  
          nodes_ex[node_id, 1] <- node_id  
          nodes_ex[node_id, 2] <- node_name  
          nodes_ex[node_id, 3] <- node1_ex
          }
          
          #matrix
          colnames(int_matrix)[node_id] <- node_name
          rownames(int_matrix)[node_id] <- node_name
          
          diet_node1 <- node1$diet_descr
          nr_food_items <- length(diet_node1)
          
          for(a in 1:nr_food_items){
            item1 <- diet_node1[[a]]
            id_item1 <- as.numeric(item1$prey_seq)  
            proportion_item1 <- as.numeric(item1$proportion)
            detritus_item1 <- as.numeric(item1$detritus_fate)
            #send to matrix
            int_matrix[id_item1,node_id] <- proportion_item1  
          }
          
        }
        
        int_matrix[is.na(int_matrix)] <- 0#replacing NA with 0
        
        if(biomass == TRUE) {  
          fwlist[[i]] <- list(biomass=nodes_biomass, pb = nodes_pb, qb = nodes_qb, ee = nodes_ee, 
            gs = nodes_gs, ex = nodes_ex, trophic_relations=int_matrix) 
        } else  { fwlist[[i]] <- int_matrix }

            #STUDY NUMBER
        if(study_id==TRUE){ fwlist[[i]]$study_id = liste_mod$model.model_number[[i]] }
        
      }
      
      names(fwlist) <- model.name
      
    })#end of outer suppressWarnings


    #REFERENCES
    if(ref==TRUE){
      
      references <- as.data.frame(matrix(ncol = 4))
      names(references) <- c("FW code", "first_author", "year", "full_ref" ) 
      
      message("Fetching the references information!")
      
      for(w in 1:nrow(l3)){
        
        #Get the reference to the vector
        full_ref1 <- as.character(l3$model.reference)[w]
        references[w,4] <- full_ref1#full reference
        references[w,1] <- as.numeric(as.character(l3$model.model_number[w]))#fw code
        references[w,2] <- as.character(l3$model.author[w])#fisrt author
        references[w,3] <- regmatches(x = full_ref1,gregexpr("[0-9]+",text = full_ref1))[[1]][1]#year
        #references[w,3] <- gsub('.+\\(([0-9]+)\\).+?$', '\\1', full_ref1)#year
        fwlist$references = references  
      }#end loop to add refs
      
      

    }#end of eb refs
    
    #ECOSYSTEM
    if(ecosyst==TRUE){
      
      ecosystem <- data.frame(l3$model.model_number, l3$model.country, l3$model.ecosystem_type)
      
      names(ecosystem) <- c("Food web", "Location", "Ecosystem")
      fwlist$ecosystem = ecosystem  

    }#end of eb ecosystem
    
    return (fwlist)
}#end of eb


#=============================================================================
#get_eb_fw 
#
# This function is meant to process a set of CSV files from a series of 
# freshwater food webs published in Christensen and Pauly, 1993. Each Lake
# has two files: Biomass... and DC... The DC file contains the dietary matrix,
# and the Biomass file contains biomass, PB, QB, and EE. Here they are 
# processed into a data structure that matches that produced by get_eb above. 
#=============================================================================

get_eb_fw = function (biomass_file_list = NULL, DC_file_list = NULL, 
    biomass = TRUE, pb = TRUE, qb = TRUE, ee = TRUE, gs = FALSE, export= FALSE, 
    study_id = FALSE, ecosyst = FALSE, ref = FALSE ) {

  fwlist = list()
  #Get model names
  model.names =  unlist(strsplit(unlist(strsplit(biomass_file_list, c("Biomass"))), c(".csv") ))
  input_list = list()
  
   for (i in 1:length(model.names) ){
        print(i)
        m2b = read.csv(paste("./fwfw/",biomass_file_list[[i]],sep=""),stringsAsFactors =F ) #get the biomass list
        colnames(m2b)[grepl("roup", names(m2b))] = "Groups"  #Make sure this category has the same name
        if(i == 12){colnames(m2b)[1] = "Groups"}
        m2b[is.na(m2b)]=0

        m2dc = read.csv(paste("./fwfw/",DC_file_list[[i]],sep=""),stringsAsFactors =F) #get the DC
        rownames(m2dc) = m2dc[,1]
        m2dc = m2dc[,-1]
        m2dc[is.na(m2dc)]=0
        m2dc = as.matrix(m2dc)

        #How many species? 
        node_names = m2b$Groups
        get_d = m2b$Groups[grepl("etritus", m2b$Groups) | grepl("Detr", m2b$Groups) ]
        get_d = get_d[!grepl("ivores", get_d)]
        if(length(get_d)<1){  #Add Detritus if it's not there
          node_names = c(node_names, "Detritus")
          d_true = rownames(m2dc)[grepl("etritus", rownames(m2dc))]
          
            nm2b = colnames(m2b)
            if("Detritus" %in% d_true == TRUE){
              d_use  = data.frame( Groups = "Detritus", biomass=sum(m2dc[grepl("etritus", rownames(m2dc)),]),
                pb= 1, qb=0)  
            }else{ 
              d_use  = data.frame( Groups = "Detritus", biomass=0,
              pb= 1, qb=0)  
            }
            m2b = full_join(m2b,d_use)
            m2b [is.na(m2b)] = 0
        }
        nnodes = dim(m2b)[1]
        nnodes2 = length(node_names)

        #Biomass, pb,qb
        nodes_biomass <- as.data.frame(matrix(ncol=3, nrow=nnodes))
        names(nodes_biomass) <- c("id", "name", "biomass")

        nodes_pb <- as.data.frame(matrix(ncol=3, nrow=nnodes))
        names(nodes_pb) <- c("id", "name", "pb")

        #qb is not always present in these data sets but we need it regardless:
        qb = "qb" %in% colnames(m2b) 
        nodes_qb <- as.data.frame(matrix(ncol=3, nrow=nnodes))
        names(nodes_qb) <- c("id", "name", "qb")
  
        #EE is not always present in these data sets but we need it regardless:
        ee = "ee"  %in% colnames(m2b) 
        nodes_ee <- as.data.frame(matrix(ncol=3, nrow=nnodes))
        names(nodes_ee) <- c("id", "name", "ee")

        for(j in 1:nnodes){
          
          node1 <- m2b[j,]
          node_id <- j 
          node_name <- node_names[j]
          
          #biomass
          node1_biomass <- as.numeric(node1$biomass)
          nodes_biomass[node_id, 1] <- node_id  
          nodes_biomass[node_id, 2] <- node_name  
          nodes_biomass[node_id, 3] <- node1_biomass  
      
          #pb
          node1_pb <- as.numeric(node1$pb)  
          nodes_pb[node_id, 1] <- node_id  
          nodes_pb[node_id, 2] <- node_name  
          nodes_pb[node_id, 3] <- node1_pb
          
          #qb
         if ( sum(qb) >0 )
          {  
          node1_qb <- as.numeric(node1$qb)  
          nodes_qb[node_id, 1] <- node_id  
          nodes_qb[node_id, 2] <- node_name  
          nodes_qb[node_id, 3] <- node1_qb 
          } else {
          #Hopefully this info is in the DC, so sum along the
          #appropriate row/column. 
          #Which is Smaller dimension? This is the one to sum across. 
          if(dim(m2dc)[1]>dim(m2dc)[2]){
             node1_qb = sum(m2dc[j,])
          } else {
             node1_qb = sum(m2dc[,j])
          }
          nodes_qb[node_id, 1] <- node_id  
          nodes_qb[node_id, 2] <- node_name  
          nodes_qb[node_id, 3] <- node1_qb 
        }
          
          if ( sum(ee) > 0 )
          {  
            node1_ee <- as.numeric(node1$ee)  
            nodes_ee[node_id, 1] <- node_id  
            nodes_ee[node_id, 2] <- node_name  
            nodes_ee[node_id, 3] <- node1_ee 
          } else {
            node1_ee = runif(1) #Uniform random number? 
            nodes_ee[node_id, 1] <- node_id  
            nodes_ee[node_id, 2] <- node_name  
            nodes_ee[node_id, 3] <- node1_ee 
          }

          
        }

        #This section is to parse out the DC. 
        #Get the orientation right. Consumers should be columns. In some cases 
        #the opposite is true and there are more columns than rows (because 
        #primary producers are not given entries) 
        if(dim(m2dc)[2]>dim(m2dc)[1]){
          m2dc=t(m2dc)
        } 

        #Normalize columns to sum to 1.
        m2dc = m2dc/matrix(colSums(m2dc),dim(m2dc)[1],dim(m2dc)[2],byrow=T )
        m2dc[!is.finite(m2dc)] = 0

        #Interaction matrix
        int_matrix = matrix(0, ncol=nnodes2, nrow=nnodes2)
        #In case dimensions are not square:
        int_matrix[1:dim(m2dc)[1], 1:dim(m2dc)[2] ] =  as.matrix(m2dc)
        colnames(int_matrix) = node_names; rownames(int_matrix) = node_names;
        int_matrix=as.data.frame(int_matrix)

        fwlist[[i]] <- list(biomass=nodes_biomass, pb = nodes_pb, qb = nodes_qb, ee = nodes_ee, 
        trophic_relations=int_matrix) 
  
        
      }
      
  names(fwlist) = model.names
  return (fwlist)

}#end of get_eb_fw




#MANGAL
get_mn = function (type = c("herbivory","predation"), ecosyst = FALSE, ref = FALSE ) {

  
  message("####################### MANGAL DATABASE #######################\n\n")
  
  message("Fetching datasets from the Mangal website! \n\n Types 'predation' and 'herbivory' by default... \n but run mangal function 'avail_type' to check available types...\n\nThis operation might take a long time!")
  
  ntypes <- length(type)
  
  net_info <- list()
  
  for(i in 1:ntypes){
    
    message(paste0("\n\nFetching information from interactions of the type ","'",type[i], "'!"))
    
    fwlist1 <- search_interactions(type = type[i]) %>% get_collection()
      
    net_info <- rbind(net_info, fwlist1)

    fwlist2 <- as.igraph(fwlist1)
    
    fwlist <- c(fwlist, fwlist2)
    
    #class(fwlist)
    
  }

  #Only keep networks that have more than presence/absence data:
  no_pa = NULL
  nnets = length(net_info)
  for(i in 1:nnets) {
    pa_tmp = net_info[[i]]$interactions$attribute.name != "presence/absence"
    no_pa = c(no_pa, sum(pa_tmp)) 
  }

  net_info = net_info[no_pa>0]
  fwlist = fwlist[no_pa] 
    
  #Converting igraph objects to data frame
  for(i in 1:length(fwlist)){
    fw2 <- fwlist[[i]]
    #convert each igraph to a data frame
    fw3 <- as_data_frame(fw2, what = "both")
    id_name <- fw3$vertices[,1:2]
    
    for(j in 1:nrow(id_name)){#clean the names
      
      node_name <- id_name$original_name[j]
      
      if (grepl(":", node_name, fixed=TRUE)) {
        node_name <- tail(strsplit(node_name, ": "))[[1]]
        id_name[j,2] <- node_name[2]
      } else id_name[j,2] <- node_name
      
      
    }#end clean names
    
    id_edges <- fw3$edges[,1:3]
    int_matrix <- as.data.frame(matrix(ncol = nrow(id_name), nrow = nrow(id_name)))
    colnames(int_matrix) <- id_name$original_name
    rownames(int_matrix) <- id_name$original_name
    
    #Fill the matrix
    for(a in 1:nrow(id_edges)){
      edge1 <- as.numeric(id_edges[a,1:2])
      name1 <- id_name[as.character(edge1[1]),][,2]
      name2 <- id_name[as.character(edge1[2]),][,2]
      int_matrix[name1,name2] <- 1
    }
    
    int_matrix[is.na(int_matrix)] <- 0 #convert all NA to zero
    
    fwlist[[i]] <- int_matrix
    
  }#end of loop to convert to a data frame

  if(ref==TRUE){
    references <- as.data.frame(matrix(ncol = 4))
    names(references) <- c("Dataset ID", "first_author", "year", "DOI" )
    
    message("Fetching references!")
    
    for(j in 1:length(net_info)){
    dataset_id <- net_info[[j]]$dataset$dataset_id
    first_author <- net_info[[j]]$reference$first_author
    year_mng <- as.numeric(net_info[[j]]$reference$year)  
    doi_mng <- net_info[[j]]$reference$doi
    references[j,1] <- dataset_id
    references[j,2] <- first_author
    references[j,3] <- year_mng
    references[j,4] <- doi_mng
    
    references <- references[order(references$`Dataset ID`),]
    rownames(references) <- 1:nrow(references)
    }
    
  }#End of mg refs
    
  if(spatial==TRUE){
    spatial1 <- as.data.frame(matrix(ncol = 4)) 
    names(spatial1) <- c("Dataset ID", "first_author", "lat", "long")
    message("Fetching coordinates!")
    
    for(z in 1: length(net_info)){
      dataset_id <- net_info[[z]]$dataset$dataset_id
      lat_mng <- net_info[[z]]$network$geom_lat
      long_mng <-  net_info[[z]]$network$geom_lon
      first_author <- net_info[[z]]$reference$first_author
      if(length(unlist(lat_mng))>1){
        
        spatial2 <- as.data.frame(matrix(ncol = 4)) 
        names(spatial2) <- c("Dataset ID", "first_author", "long", "lat" )
        
        for(b in 1:length(unlist(lat_mng))){
          spatial2[b,3] <- long_mng[[1]] [b]
          spatial2[b,4] <- lat_mng [[1]] [b]
          }
        
        spatial2[,1] <- dataset_id
        spatial2[,2] <- first_author
        
        spatial1 <- rbind(spatial1, spatial2)
        
      }
      spatial1[z,1] <- dataset_id
      spatial1[z,2] <- first_author
      if(length(unlist(lat_mng))==1) spatial1[z,3] <- lat_mng
      if(length(unlist(lat_mng))==1) spatial1[z,4] <- long_mng
    }
      
      spatial1 <- spatial1[order(spatial1$`Dataset ID`),]
      rownames(spatial1) <- 1:nrow(spatial1)
      
  }#End of mg spatial
    
    if (exists("references") & exists("spatial1")) (if(nrow(references)!=nrow(spatial1)) message("WARNING: There are more than on FW in some datasets! References and Spatial data frames have different number of rows."))

return (fwlist)
    
}#end of mangal