# Differential Expression Analysis 
# COSMOS pipeline. Goal: Create mechanistic networks based on CLL and HD LOG2FC values. 

library(readxl)
library(readr)
library(biomaRt)
library(dplyr)
library(reshape2)
library(pheatmap)
library(vsn)
library(stringr)

############ Paths to data AMC PC ############ 
### Getting DGE for RNA
## Gene expression
# T0
T0_DEG_path <- "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\Transcriptomic\\EC024_output_data_T0.RData"
load(T0_DEG_path)
write.csv(DEG_df,file = "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\Transcriptomic\\T0_DEG.csv")

# T48
rm(list = ls())
T48_DEG_path <- "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\Transcriptomic\\EC024_output_data_T48.RData"

load(T48_DEG_path)
write.csv(DEG_df,file = "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\Transcriptomic\\T48_DEG.csv")


###
rm(list = ls())

########################################## paths frfr for both RNA and metabs

T0_DEG_path <- "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\Transcriptomic\\T0_DEG.csv"
T48_DEG_path <- "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\Transcriptomic\\T48_DEG.csv"

## Metabolites
T0_metab_path <- "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\Metabolomic\\T0_Statistical_Data.xlsx"
T48_metab_path <- "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\Metabolomic\\T48_Statistical_Data.xlsx"

metabolite_matching <- as.data.frame(read_excel("L:/basic/divg/EXIM/ImmunoHematology/Cathy Magnée/Data/CD3_HDvsCLL/Metabolomic/Statistical Data_Incl_HMDB.xlsx"))


########################################## Data exploration - RNA
T0_DEG <- as.data.frame(read_csv(T0_DEG_path))

T0_DEG <- T0_DEG[, (colnames(T0_DEG) %in% c("GENE", "log2FoldChange"))]

T48_DEG <- as.data.frame(read_csv(T48_DEG_path))

T48_DEG <- T48_DEG[, (colnames(T48_DEG) %in% c("GENE", "log2FoldChange"))]

hist(as.numeric(unlist(T0_DEG$log2FoldChange)), breaks = 18130)
hist(as.numeric(unlist(T48_DEG$log2FoldChange)), breaks = 17057)

hist(as.numeric(unlist(T0_DEG$log2FoldChange)), breaks = 18130, xlim = c(-3,3))
hist(as.numeric(unlist(T48_DEG$log2FoldChange)), breaks = 17057, xlim = c(-3,3))


########################################## Data exploration - metab
T0_DA_metab <- as.data.frame(read_excel(T0_metab_path))
T48_DA_metab <- as.data.frame(read_excel(T48_metab_path))

T0_DA_metab <- T0_DA_metab[, (colnames(T0_DA_metab) %in% c("HMDB.id", "Log2_fold_change"))]
T48_DA_metab <- T48_DA_metab[, (colnames(T48_DA_metab) %in% c("HMDB.id", "Log2_fold_change"))]

hist(as.numeric(unlist(T0_DA_metab$Log2_fold_change)), breaks = 119)
hist(as.numeric(unlist(T48_DA_metab$Log2_fold_change)), breaks = 121)

########################### Explanation
# Now we have, for both T0 and T48, the differential analysis for gene expression and metabolites 
# Lets infer tf activity from this?

# ____________________________________________________________________________________________
# 4 - Preparing COSMOS inputs
# ____________________________________________________________________________________________


########## Getting collectri database 
library(readr)
library(decoupleR)
library(OmnipathR)
library(dplyr)

load("L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\CosmosR\\COSMOS_Collectri_db.RData")
#################### THIS IS HOW I MAKE THE ABOVE DF
# # Omnipath settings
# collectri_regulons <- get_collectri(organism = 'human', split_complexes = FALSE)
# collectri_regulons <- decoupleR::get_collectri()
# collectri_regulons <- collectri_regulons[!(is.na(collectri_regulons[, 2]) | collectri_regulons[, 2] == ""), ]
# # Omnipath settings to get collectri regulons
# 
# omnipath_set_cachedir('tmpcache')
# omnipath_set_console_loglevel('trace')
# 
# # Getting all regulons
# ci <- collectri() # overview
# 
# ci_regulons <- ci[,c(3,4,6)] # only need these columns from overview. These are source, target and stimulates (0 and 1)
# 
# colnames(ci_regulons) <- c("source", "target", "mor")
# ci_regulons <- unique(ci_regulons) # remove pure dupes
# 
# # Some of the regulons are duplicated with different mor.
# # This means the same source/target combination has both 0 and 1 options
# # Let's change this.
# temp <- ci[duplicated(ci[, 3:4]) | duplicated(ci[, 3:4], fromLast = TRUE), ] # all duplicates, also with the same mor values
# duplicates <- anti_join(temp, ci[duplicated(ci[, c(3,4,6)]) | duplicated(ci[, c(3,4,6)], fromLast = TRUE), ]) # only different mor values
# duplicates <- duplicates[,c(3,4,6)]
# colnames(duplicates) <- c("source", "target", "mor")
# duplicates <- unique(duplicates)
# 
# # For loop that checks the value of the duplicates in ci_regulons what their values are in collectri_regulons then changes accordingly
# for(i in 1:nrow(duplicates)){
#   source <- duplicates$source[i]
#   target <- duplicates$target[i]
# 
#   collectri_row <- collectri_regulons %>%
#     filter(source == !!source & target == !!target)
# 
#   if (nrow(collectri_row) > 0){ # if it is in collectri
#     # check mor value
#     if (collectri_row$mor == -1){
#       ci_regulons <- ci_regulons |>
#         mutate(mor = ifelse(source == !!source & target == !!target, 0 , mor))
#     } else {
#       ci_regulons <- ci_regulons |>
#         mutate(mor = ifelse(source == !!source & target == !!target, 1, mor))
#     }
#   }
# }
# 
# ci_regulons <- unique(ci_regulons)
# 
# 
# # Change inhibition relation to -1
# ci_regulons$mor[ci_regulons$mor == 0] <- -1
# 
# save(ci_regulons, file = "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\CosmosR\\COSMOS_Collectri_db.RData")

########################## ULM ##########################
T0_DEG <- unique(T0_DEG)
T0_DEG <- T0_DEG[!(duplicated(T0_DEG$GENE)), ]
T0_DEG <- T0_DEG[!(is.na(T0_DEG$GENE)), ]

RNA_scaled <- T0_DEG
rownames(RNA_scaled) <- RNA_scaled$GENE
RNA_scaled$GENE <- NULL

T0_TF_activities <- apply(RNA_scaled,2,function(x){
  x <- as.data.frame(x[which(!is.na(x))]) # no NAs allowed
  TFs <- run_ulm(as.matrix(x), # actual running of TF activity estimation function
                 network = unique(ci_regulons),
                 minsize = 20)
  TFs <- as.data.frame(TFs) # Reformatting into cosmos input format
  TFs <- TFs[which(TFs$statistic == "ulm"),c(2,4)]
  as_input <- TFs[,2]
  names(as_input) <- TFs[,1]
  return(as_input)
})


T48_DEG <- unique(T48_DEG)
T48_DEG <- T48_DEG[!(duplicated(T48_DEG$GENE)), ]
T48_DEG <- T48_DEG[!(is.na(T48_DEG$GENE)), ]

RNA_scaled <- T48_DEG
rownames(RNA_scaled) <- RNA_scaled$GENE
RNA_scaled$GENE <- NULL

T48_TF_activities <- apply(RNA_scaled,2,function(x){
  x <- as.data.frame(x[which(!is.na(x))]) # no NAs allowed
  TFs <- run_ulm(as.matrix(x), # actual running of TF activity estimation function
                 network = unique(ci_regulons),
                 minsize = 20)
  TFs <- as.data.frame(TFs) # Reformatting into cosmos input format
  TFs <- TFs[which(TFs$statistic == "ulm"),c(2,4)]
  as_input <- TFs[,2]
  names(as_input) <- TFs[,1]
  return(as_input)
})  

# ____________________________________________________________________________________________
# Prep data for Cosmos_Inputs
# ____________________________________________________________________________________________

T0_transcriptomic <- T0_DEG[!(is.na(T0_DEG$GENE)), ]
row.names(T0_transcriptomic) <- T0_transcriptomic$GENE
T0_transcriptomic$GENE <- NULL
T0_transcriptomic <- apply(T0_transcriptomic, 2, function(x){
  return(x[which(!is.na(x))])
}, simplify = F)

T48_transcriptomic <- T48_DEG[!(is.na(T48_DEG$GENE)), ]
row.names(T48_transcriptomic) <- T48_transcriptomic$GENE
T48_transcriptomic$GENE <- NULL
T48_transcriptomic <- apply(T48_transcriptomic, 2, function(x){
  return(x[which(!is.na(x))])
}, simplify = F)


T0_metabolomic <- T0_DA_metab[!(is.na(T0_DA_metab$HMDB.id)), ]
row.names(T0_metabolomic) <- T0_metabolomic$HMDB.id
T0_metabolomic$HMDB.id <- NULL
colnames(T0_metabolomic) <- "log2FoldChange"
T0_metabolomic <- apply(T0_metabolomic, 2, function(x){
  x <- x[which(!is.na(x))]
  return(x)
},simplify = F)

T48_metabolomic <- T48_DA_metab[!(is.na(T48_DA_metab$HMDB.id)), ]
row.names(T48_metabolomic) <- T48_metabolomic$HMDB.id
T48_metabolomic$HMDB.id <- NULL
colnames(T48_metabolomic) <- "log2FoldChange"
T48_metabolomic <- apply(T48_metabolomic, 2, function(x){
  x <- x[which(!is.na(x))]
  return(x)
},simplify = F)

T0_TF_activities <- apply(T0_TF_activities, 2, function(x){
  x <- x[which(!is.na(x))]
  return(x)
},simplify = F)
T48_TF_activities <- apply(T48_TF_activities, 2, function(x){
  x <- x[which(!is.na(x))]
  return(x)
},simplify = F)
# ____________________________________________________________________________________________
# Create Cosmos_Inputs
# ____________________________________________________________________________________________

TF_activities <- T0_TF_activities
RNA_scaled <- T0_transcriptomic
metabs_scaled <- T0_metabolomic


T0_cosmos_inputs <- sapply(names(TF_activities), function(cell_line, RNA_scaled, metabs_scaled, TF_activities){
  cosmos_input <- list()
  cosmos_input[["RNA"]] <- RNA_scaled[[cell_line]]
  cosmos_input[["metabolomic"]] <- metabs_scaled[[cell_line]]
  cosmos_input[["TF_scores"]] <- TF_activities[[cell_line]]
  return(cosmos_input)
}, RNA_scaled = RNA_scaled, metabs_scaled = metabs_scaled, TF_activities= TF_activities, USE.NAMES = T, simplify = F)

TF_activities <- T48_TF_activities
RNA_scaled <- T48_transcriptomic
metabs_scaled <- T48_metabolomic


T48_cosmos_inputs <- sapply(names(TF_activities), function(cell_line, RNA_scaled, metabs_scaled, TF_activities){
  cosmos_input <- list()
  cosmos_input[["RNA"]] <- RNA_scaled[[cell_line]]
  cosmos_input[["metabolomic"]] <- metabs_scaled[[cell_line]]
  cosmos_input[["TF_scores"]] <- TF_activities[[cell_line]]
  return(cosmos_input)
}, RNA_scaled = RNA_scaled, metabs_scaled = metabs_scaled, TF_activities= TF_activities, USE.NAMES = T, simplify = F)

save(T0_cosmos_inputs, file = "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\CosmosR\\T0_cosmos_inputs_DA.RData")
save(T48_cosmos_inputs, file = "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\CosmosR\\T48_cosmos_inputs_DA.RData")


# ____________________________________________________________________________________________
# Run COSMOS
# ____________________________________________________________________________________________

# Meta stuff
data("meta_network")

meta_network <- meta_network_cleanup(meta_network)
my_options <- default_CARNIVAL_options(solver = "cplex")
my_options$solverPath <- "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64"
my_options$workdir <- "C:/Users/cmagnee/Documents/cplex"
my_options$outputFolder <- "C:/Users/cmagnee/Documents/cplex"
my_options$solver <- "cplex"
my_options$mipGAP <- 0.05
my_options$poolCap <- 500
my_options$limitPop <- 2500
my_options$threads <- 19 # Adjust based on your machine
# my_options$threads <- 6 # My laptop

max_depth <- 3
seconds_per_step <- 1 * 3600 / 64 # Time in h * Seconds in an hour / Amount of steps 
my_options$timelimit <- seconds_per_step  
time_df <- data.frame()


inputs <- data.frame()
# cosmos_inputs <- T0_cosmos_inputs
cosmos_inputs <- T48_cosmos_inputs # T0 or T48
for (patient in names(cosmos_inputs)){
  
  sig_input <- cosmos_inputs[[patient]]$TF_scores
  metab_input <- cosmos_inputs[[patient]]$metabolomic

  RNA_input <- cosmos_inputs[[patient]]$RNA
  metab_input <- prepare_metab_inputs(metab_input, c("c","m"))

  sig_input <- sig_input[abs(sig_input) > 3.5]
  metab_input <- metab_input[abs(metab_input) > 0.7]

  metab_input <- cosmosR:::filter_input_nodes_not_in_pkn(metab_input, meta_network)
  sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)
  new_row <- list(
    patient = patient,
    n_sig_inpul = length(sig_input),
    n_metab_input = length(metab_input)
  )
  inputs <- rbind(inputs, new_row)
}
  
start.time <- Sys.time()
for (i in 2:2){
  cosmos_inputs <- list(T0_cosmos_inputs, T48_cosmos_inputs)[[i]]
  stim_point <- list("T0", "T48")[[i]]
  for (patient in names(cosmos_inputs)){
    skip_to_next <<- F
    tryCatch({
      sig_input <- cosmos_inputs[[patient]]$TF_scores
      metab_input <- cosmos_inputs[[patient]]$metabolomic 
      RNA_input <- cosmos_inputs[[patient]]$RNA
      
      # Filter inputs
      metab_input <- prepare_metab_inputs(metab_input, c("c","m"))
      
      
      filter <- "s2m2" ## This one is important!! s = signalling input, m = metabolite input filter do not put in the points or smtng
      sig_input <- sig_input[abs(sig_input) > 3.5] 
      metab_input <- metab_input[abs(metab_input) > 0.7] 
      
      metab_input <- cosmosR:::filter_input_nodes_not_in_pkn(metab_input, meta_network)
      sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)
      
      # Stp 1
      test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = meta_network[,c(1,3,2)],
                                                            signaling_data = sig_input,
                                                            metabolic_data = metab_input,
                                                            diff_expression_data = RNA_input,
                                                            maximum_network_depth = max_depth,
                                                            remove_unexpressed_nodes = T,
                                                            filter_tf_gene_interaction_by_optimization = T,
                                                            CARNIVAL_options = my_options)
      my_options$timelimit <- seconds_per_step
      new_row <- list(
        patient = patient,
        step = "step_1",
        end_time = round(Sys.time()-start.time,2)
      )
      time_df <- rbind(time_df, new_row)
      #Stp 2
      test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                            CARNIVAL_options = my_options)
      new_row <- list(
        patient = patient,
        step = "step_2",
        end_time = round(Sys.time()-start.time,2)
      )
      time_df <- rbind(time_df, new_row)
      
      formatted_res <- format_COSMOS_res(test_result_for)
      
      SIF <- formatted_res[[1]]
      ATT <- formatted_res[[2]]
      
      SIF <- SIF[which(SIF$Weight != 0),]
      
      RNA_input_df <- data.frame(Nodes = names(RNA_input), t = RNA_input)
      ATT <- merge(ATT, RNA_input_df, all.x = T)
      ATT <- ATT[ATT$AvgAct != 0,]
      
      
      # write_csv(SIF, file = paste("../results/",paste(patient, "_SIF.csv",sep = ""), sep = ""))
      # write_csv(ATT, file = paste("../results/",paste(patient, "_ATT.csv",sep = ""), sep = ""))
      
      #Stp 3
      my_options$timelimit <- seconds_per_step
      
      
      test_back <- preprocess_COSMOS_metabolism_to_signaling(meta_network = meta_network[,c(1,3,2)],
                                                             signaling_data = sig_input,
                                                             metabolic_data = metab_input,
                                                             diff_expression_data = RNA_input,
                                                             maximum_network_depth = max_depth,
                                                             remove_unexpressed_nodes = TRUE,
                                                             filter_tf_gene_interaction_by_optimization = TRUE,
                                                             CARNIVAL_options = my_options)
      new_row <- list(
        patient = patient,
        step = "step_3",
        end_time = round(Sys.time()-start.time,2)
      )
      time_df <- rbind(time_df, new_row)
      #Stp 4
      my_options$timelimit <- seconds_per_step
      
      test_result_back <- run_COSMOS_metabolism_to_signaling(data = test_back,
                                                             CARNIVAL_options = my_options)
      new_row <- list(
        patient = patient,
        step = "step_4",
        end_time = round(Sys.time()-start.time,2)
      )
      time_df <- rbind(time_df, new_row)
      formatted_res_back <- format_COSMOS_res(test_result_back)
      
      SIF_back <- formatted_res_back[[1]]
      ATT_back <- formatted_res_back[[2]]
      
      SIF_back <- SIF_back[which(SIF_back$Weight != 0),]
      
      ATT_back <- merge(ATT_back, RNA_input_df, all.x = T)
      ATT_back <- ATT_back[ATT_back$AvgAct != 0,]
      
      # write_csv(SIF_back, file = paste("../results/",paste(patient, "_SIF_back.csv",sep = ""), sep = ""))
      # write_csv(ATT_back, file = paste("../results/",paste(patient, "_ATT_back.csv",sep = ""), sep = ""))
      
      #Stp 5: Merge all & save
      
      SIF_full <- as.data.frame(rbind(SIF,SIF_back))
      SIF_full <- unique(SIF_full)
      
      ATT_full <- as.data.frame(rbind(ATT,ATT_back))
      ATT_full <- unique(ATT_full)
      
      P_nodes <- ATT_full[ATT_full$NodeType == "P","Nodes"]
      M_nodes <- ATT_full[ATT_full$NodeType == "M","Nodes"]
      C_nodes <- union(P_nodes,M_nodes)
      # C_nodes <- C_nodes[which(C_nodes %in% c(SIF$Node1,SIF$Node2) & C_nodes %in% c(SIF_back$Node1,SIF_back$Node2))]
      
      ATT_full <- ATT_full[,-6]
      
      ATT_full <- ATT_full %>% group_by(Nodes) %>% summarise_each(funs(mean(., na.rm = TRUE)))
      ATT_full <- as.data.frame(ATT_full)
      
      ATT_full$NodeType <- ifelse(ATT_full$Nodes %in% C_nodes,"C",ifelse(ATT_full$Nodes %in% P_nodes,"P",ifelse(ATT_full$Nodes %in% M_nodes,"M","")))
      
      # Automatically create path to save results, so it does not get messy in my folders
      core_path <- "L:/basic/divg/EXIM/ImmunoHematology/Cathy Magnée/Data/CD3_HDvsCLL/CosmosR"
      
      result_path <- paste("/Networks/d", max_depth, sep = "") # Add depth folder
      dir.create(file.path(core_path, result_path), showWarnings = FALSE)
      result_path <- paste(result_path, paste("/f", filter, sep = ""), sep = "") # Add filter folder
      dir.create(file.path(core_path, result_path), showWarnings = FALSE)
      result_path <- paste(result_path, paste("/", stim_point, sep = ""), sep = "") # Add stim/unstim point folder
      dir.create(file.path(core_path, result_path), showWarnings = FALSE)
      result_path <- paste(result_path, "/DA/", sep = "") # Add csv folder
      dir.create(file.path(core_path, result_path), showWarnings = FALSE)
      
      result_path <- paste(core_path, result_path, sep = "")
      # result_path <- paste("../results/Network_files_csvs/d", paste(max_depth, paste("/f", paste (filter, paste("/"), paste(stim_point, "/", sep = ""), sep = ""), sep = ""), sep = ""), sep = "")
      write_csv(SIF_full, file = paste(result_path, paste(patient, "_SIF_full.csv", sep = ""), sep = ""))
      write_csv(ATT_full, file = paste(result_path, paste(patient, "_ATT_full.csv", sep = ""), sep = ""))
      
      # How long did this take?
      new_row <- list(
        patient = patient,
        step = 'final',
        end_time = round(Sys.time()-start.time,2))
      time_df <- rbind(time_df, new_row)
    }, error = function(e) {
      message(paste("Error in processing patient:", patient, "-", e$message))
      skip_to_next <<- TRUE  # Set flag to skip the rest of the iteration
    })
    if (skip_to_next){
      next
    }
  }
}
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken


# ____________________________________________________________________________________________
# 6 - Preprocess network
# ____________________________________________________________________________________________

for (i in 2:2){
  stim_point <- list("T0", "T48")[[i]]
  # Path to where the data should be at AMC pc:
  core_path <- "L:/basic/divg/EXIM/ImmunoHematology/Cathy Magnée/Data/CD3_HDvsCLL/CosmosR"
  
  result_path <- paste("/Networks/d", max_depth, sep = "") # Add depth folder
  result_path <- paste(result_path, paste("/f", filter, sep = ""), sep = "") # Add filter folder
  result_path <- paste(result_path, paste("/", stim_point, sep = ""), sep = "") # Add stim/unstim point folder
  result_path <- paste(result_path, "/DA/", sep = "") # Add csv folder
  
  data_path <- paste(core_path, result_path, sep = "")
  
  # Edit below to what you want to preprocess
  all_files <- list.files(path = data_path, pattern = "*_SIF_full.csv", full.names = TRUE)
  for (files in all_files){
    patient <- str_extract(files, "[^/]+(?=_SIF)")
    
    # Importing SIF and ATT file for each available patient / HD sample
    SIF_path <- paste(data_path, paste(patient, "_SIF_full.csv", sep = ""), sep = "")
    SIF <- read_csv(SIF_path)
    SIF <- unique(SIF)
    
    ATT_path <- paste(data_path, paste(patient, "_ATT_full.csv",sep = ""), sep = "")
    ATT <- read_csv(ATT_path)
    ATT <- unique(ATT)
    
    # Creating MoleType df for Metabs 
    MolType_df <- data.frame()
    for(nodes in ATT$Nodes){
      # print(substr(nodes, 0, 5))
      if(substr(nodes, 0, 5) == "Metab"){
        new_row <- list(Node = nodes,
                        MoleculeType = "Metabolite")
        MolType_df <- rbind(MolType_df, new_row)
      } else if(substr(nodes, 0, 5) == "Enzym"){
        new_row <- list(Node = nodes,
                        MoleculeType = "Enzyme")
        MolType_df <- rbind(MolType_df, new_row)
      } else {
        new_row <- list(Node = nodes,
                        MoleculeType = NA)
        MolType_df <- rbind(MolType_df, new_row)
      }
      
    }
    # Creating MoleculeType column. This will be built on in 7 - Cytoscape.Rmd
    ATT <- merge(ATT, MolType_df, by.x = "Nodes", by.y = "Node", all.x = T)
    
    # Sub Enzyme0000__ (random numbers) to Enzyme_
    SIF$Node1 <- sub("Enzyme[0-9]+__", "Enzyme_", SIF$Node1)
    SIF$Node2 <- sub("Enzyme[0-9]+__", "Enzyme_", SIF$Node2)
    ATT$Nodes <- sub("Enzyme[0-9]+__", "Enzyme_", ATT$Nodes)
    
    #Filter these rando Enzymes to collaborate (all become 1, mean of availables)
    # ATT <- ATT %>%
    #   group_by(Nodes) %>% # Group by the Nodes column
    #   summarize(across(everything(), mean, na.rm = TRUE))
    
    # Remove Metab__ part of Node names
    SIF$Node1 <- sub(".*__", "", SIF$Node1)
    SIF$Node2 <- sub(".*__", "", SIF$Node2)
    ATT$Nodes<- sub(".*__", "", ATT$Nodes)
    
    SIF$Node1 <- sub(',', "", SIF$Node1)
    SIF$Node2 <- sub(',', "", SIF$Node2)
    ATT$Nodes<- sub(',', "", ATT$Nodes)
    SIF$Node1 <- sub(',', "", SIF$Node1)
    SIF$Node2 <- sub(',', "", SIF$Node2)
    ATT$Nodes<- sub(',', "", ATT$Nodes)
    
    SIF <- SIF |>
      filter(Node1 != Node2) |> # First removing nodes that are re-activating themselves
      group_by(Node1, Sign, Node2) |>
      summarize(Weight = sum(Weight), .groups='drop')
    
    # Remove _c and _m from Node names (metabolites). But write down in NodeType # Cancelled because it gives information about cell system.
    # ATT$NodeType <- ifelse(grepl("_c$", ATT$Nodes) & is.na(ATT$NodeType), "C", ATT$NodeType)
    # ATT$NodeType <- ifelse(grepl("_m$", ATT$Nodes), "M", ATT$NodeType)
    # 
    # SIF$Node1 <- sub("_m", "", SIF$Node1)
    # SIF$Node2 <- sub("_m", "", SIF$Node2)
    # ATT$Nodes<- sub("_m", "", ATT$Nodes)
    # 
    # SIF$Node1 <- sub("_c", "", SIF$Node1)
    # SIF$Node2 <- sub("_c", "", SIF$Node2)
    # ATT$Nodes<- sub("_c", "", ATT$Nodes)
    
    # # Add additional column "Acticity_ (patient network name)". This will make network processing within cytoscape easier.
    # temp <- paste("Activity_", patient, "")
    # ATT[temp] <- ATT$Activity
    # 
    # # Add additional column "In_(patient network name)". This will make network processing within cytoscape easier.
    # temp <- paste("In_", patient,sep = "")
    # ATT[temp] <- 1
    
    # Saving processed csvs
    # result_path <- paste("/processed_results/d", max_depth, sep = "") # Laptop only
    data_path <- paste(core_path, result_path, sep = "") # Laptop only
    write_csv(SIF, file = paste(data_path, paste(patient, "_SIF_processed.csv", sep = ""), sep = ""))
    write_csv(ATT, file = paste(data_path, paste(patient, "_ATT_processed.csv", sep = ""), sep = ""))
    
  }
}



# ____________________________________________________________________________________________
# 7 - Cytoscape
# ____________________________________________________________________________________________

style.name = "empty"
Create_Network <- function(ATT, SIF, MergedGroup){
  # Translate this to Node Table and Edges Table for in Cytoscape
  nodes <- data.frame(id = ATT$Nodes,
                      Activity = ATT$Activity,
                      measured = ATT$measured,
                      moleculeType = ATT$MoleculeType,
                      stringsAsFactors = F)
  edges <- data.frame(source = SIF$Node1,
                      target = SIF$Node2,
                      interaction = as.character(SIF$Sign),
                      weight = SIF$Weight,
                      stringsAsFactors = FALSE)
  
  # Two things: 1: Adding moleculeType column to nodes table; 2: putting shapes to these
  # Step 1 should be done earlier, maybe in 6 - Preprocess_network.Rmd.
  #StringDB Tryout
  ATT_string_interaction_cmd <- paste('string protein query taxonID=9606 cutoff=0.99 query=', paste(nodes$id, collapse=","),'"',sep="")
  commandsGET(ATT_string_interaction_cmd )
  string_nodes <- getTableColumns(table = "node", columns = c("display name", "target::family"), network = "STRING network")
  colnames(string_nodes) <- c("id", "moleculeType")
  write.csv(string_nodes, file = paste(data_path, paste(MergedGroup, ".csv", sep = ""), sep = "StringDB_"))
  string_nodes <- merge(nodes, string_nodes, by = "id", all.x = T)
  nodes <- string_nodes |>
    mutate(moleculeType = coalesce(moleculeType.x, moleculeType.y, "Other"))
  nodes$moleculeType.x <- NULL
  nodes$moleculeType.y <- NULL
  
  # closeSession(save.before.closing = F)
  
  
  createNetworkFromDataFrames(nodes, edges, title = MergedGroup, collection = MergedGroup) # if(substr(MergedGroup, 0, 3)== "CLL"){"CLL"} else { "HD"})
  # You should see a network rn frfr
  
  # Setting style for network OR creating it first, then setting the style. If else so we do not remake this style over and over.
  if (style.name == "COSMOS_Style"){
    setVisualStyle(style.name)
  } else {
    style.name = "COSMOS_Style"
    defaults <- list(NODE_SHAPE="ROUND_RECTANGLE",
                     NODE_SIZE=50,
                     EDGE_TRANSPARENCY=255)
    nodeLabels <- mapVisualProperty('node label','id','p')
    nodeFills <- mapVisualProperty('node fill color',
                                   'Activity',
                                   'c',
                                   c(-1, 1),
                                   c('#CC0000', '#009900'))
    nodeBorderWidth <- mapVisualProperty('node border width',
                                         'measured',
                                         'd',
                                         c(0.0, 1.0),
                                         c(0, 4))
    arrowShapes <- mapVisualProperty('Edge Target Arrow Shape', # What style property we are mapping
                                     'interaction', # What table column we are using for it
                                     'd', # Type of mapping (discrete)
                                     c("-1", "1"),c("T","Arrow"))
    edgeWidth <- mapVisualProperty('edge width','weight','p')
    
    
    createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills,nodeBorderWidth,arrowShapes,edgeWidth))
    setVisualStyle(style.name)
  }
  
  
  # setLayoutProperties("degree-circle")
  # 2: Connecting shapes to moleculeType
  
  column <- 'moleculeType'
  values <- c ('Metabolite', 'Kinase',
               'Transcription Factor', 'Other',
               'Enzyme', 'Epigenetic',
               'GPCR', 'TF-Epigenetic',
               'Nuclear Receptor', 'Transporter')
  shapes <- c ('ELLIPSE', 'TRIANGLE',
               'DIAMOND', 'PARALLELOGRAM',
               'RECTANGLE', 'VEE',
               'HEXAGON', 'DIAMOND',
               'OCTAGON', 'ROUND RECTANGLE')
  setNodeShapeMapping (column, values, shapes, style.name = style.name)
  layoutNetwork("degree-circle")
  
  # Save session & a pdf picture of it
  saveSession(filename = paste(data_path, MergedGroup, sep = "Merged_Network_"))
  
  exportImage(paste(data_path, MergedGroup, sep = "Degree_Circle_"), 'PDF') #.pdf
  
  
  # Select nodes that have a high interconnectivity (>10 edges to other nodes)
  createDegreeFilter("ImportantNodes", criterion = c(10,100), network = MergedGroup, apply = T)
  Important_Nodes <- data.frame(getSelectedNodes())
  Important_Nodes <- na.exclude(Important_Nodes)
  
  # IF there are highly interconnected nodes, safe them and their activity.
  
}

for (i in 1:2){
  cosmos_inputs <- list(T0_cosmos_inputs, T48_cosmos_inputs)[[i]]
  stim_point <- list("T0", "T48")[[i]]
  # Path to where the data should be at AMC pc:
  core_path <- "L:/basic/divg/EXIM/ImmunoHematology/Cathy Magnée/Data/CD3_HDvsCLL/CosmosR"
  
  result_path <- paste("/Networks/d", max_depth, sep = "") # Add depth folder
  result_path <- paste(result_path, paste("/f", filter, sep = ""), sep = "") # Add filter folder
  result_path <- paste(result_path, paste("/", stim_point, sep = ""), sep = "") # Add stim/unstim point folder
  result_path <- paste(result_path, "/DA/", sep = "") # Add csv folder
  
  data_path <- paste(core_path, result_path, sep = "")
  for (patient in names(cosmos_inputs)){
    SIF_path <- paste(data_path, paste(patient, "_SIF_processed.csv", sep = ""), sep = "")
    SIF <- read_csv(SIF_path, show_col_types = F)
    ATT_path <- paste(data_path, paste(patient, "_ATT_processed.csv",sep = ""), sep = "")
    ATT <- read_csv(ATT_path, show_col_types = F)
    # 
    Create_Network(ATT = ATT, SIF= SIF, MergedGroup = patient)
  }
}
# cosmos_inputs <- T48_cosmos_inputs
# stim_point <- "T48"
# core_path <- "L:/basic/divg/EXIM/ImmunoHematology/Cathy Magnée/Data/CD3_HDvsCLL/CosmosR"
# 
# result_path <- paste("/Networks/d", max_depth, sep = "") # Add depth folder
# result_path <- paste(result_path, paste("/f", filter, sep = ""), sep = "") # Add filter folder
# result_path <- paste(result_path, paste("/", stim_point, sep = ""), sep = "") # Add stim/unstim point folder
# result_path <- paste(result_path, "/DA/", sep = "") # Add csv folder
# 
# data_path <- paste(core_path, result_path, sep = "")
# 
# SIF_path <- paste(data_path, paste(patient, "_SIF_processed.csv", sep = ""), sep = "")
# SIF <- read_csv(SIF_path, show_col_types = F)
# ATT_path <- paste(data_path, paste(patient, "_ATT_processed.csv",sep = ""), sep = "")
# ATT <- read_csv(ATT_path, show_col_types = F)
# # 
# 
# nodes <- data.frame(id = ATT$Nodes,
#                     Activity = ATT$Activity,
#                     measured = ATT$measured,
#                     moleculeType = ATT$MoleculeType,
#                     stringsAsFactors = F)
# edges <- data.frame(source = SIF$Node1,
#                     target = SIF$Node2,
#                     interaction = as.character(SIF$Sign),
#                     weight = SIF$Weight,
#                     stringsAsFactors = FALSE)
# 
# # Two things: 1: Adding moleculeType column to nodes table; 2: putting shapes to these
# # Step 1 should be done earlier, maybe in 6 - Preprocess_network.Rmd.
# #StringDB Tryout
# ATT_string_interaction_cmd <- paste('string protein query taxonID=9606 cutoff=0.99 query=', paste(nodes$id, collapse=","),'"',sep="")
# commandsGET(ATT_string_interaction_cmd )
# string_nodes <- getTableColumns(table = "node", columns = c("display name", "target::family"), network = "STRING network")
# colnames(string_nodes) <- c("id", "moleculeType")
# write.csv(string_nodes, file = paste(data_path, paste(MergedGroup, ".csv", sep = ""), sep = "StringDB_"))
# string_nodes <- merge(nodes, string_nodes, by = "id", all.x = T)
# nodes <- string_nodes |>
#   mutate(moleculeType = coalesce(moleculeType.x, moleculeType.y, "Other"))
# nodes$moleculeType.x <- NULL
# nodes$moleculeType.y <- NULL
# 
# 
# createNetworkFromDataFrames(nodes, edges, title = stim_point, collection = stim_point) 
