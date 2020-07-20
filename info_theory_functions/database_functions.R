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
get_eb = function (biomass = TRUE, pb = TRUE, qb = TRUE, ee = TRUE, ecosyst = FALSE, ref = FALSE ) {

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
      }#end of loop to get input list
      
      mnames <- names(input_list)
      
      for (i in 1:length(input_list)){
        
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
           names(nodes_ee) <- c("id", "name", "qb")
           }
        
        
        int_matrix <- as.data.frame(matrix(ncol=nnodes, nrow=nnodes))
        
        for(j in 1:length(m2)){
          
          node1 <- m2[[j]]
          node_id <- as.numeric(node1$group_seq)  
          node1_biomass <- as.numeric(node1$biomass)  
          node1_pb <- as.numeric(node1$pb)  
          node1_qb <- as.numeric(node1$qb)  
          node1_ee <- as.numeric(node1$ee)  

          node_name <- node_names[j]
          
          #biomass
          if (biomass == TRUE)
          {  
          nodes_biomass[node_id, 1] <- node_id  
          nodes_biomass[node_id, 2] <- node_name  
          nodes_biomass[node_id, 3] <- node1_biomass  
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
          fwlist[[i]] <- list(biomass=nodes_biomass, trophic_relations=int_matrix) 
        } else  { fwlist[[i]] <- int_matrix
        
        }
        
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
