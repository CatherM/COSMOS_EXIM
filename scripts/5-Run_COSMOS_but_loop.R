setwd("L:/basic/divg/EXIM/ImmunoHematology/Cathy Magnée/Scripts/COSMOS_EXIM/scripts")
start.time <- Sys.time()
library(cosmosR)
library(readr)
library(dplyr)

T0_input_path <- "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\CosmosR\\T0_cosmos_inputs.RData"
T48_input_path <- "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\CosmosR\\T48_cosmos_inputs.RData"

data("meta_network")

meta_network <- meta_network_cleanup(meta_network)

load(T0_input_path)
load(T48_input_path)

cosmos_inputs <- T48_cosmos_inputs # T0 or T48
names(cosmos_inputs)

max_depth <- 3
seconds_per_step <- 3600 * 1.25

#In order to adapt options to users specification we can load them into a variable 
#that will then be passed to preprocess_COSMOS_signaling_to_metabolism CARNIVAL_options parameter
my_options <- default_CARNIVAL_options(solver = "cplex")

#Here the user should provide a path to its CPLEX executable

# my_options$solverPath <- "~/Documents/cplex" #or cbc solver executable

my_options$solverPath <- "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64"

# Change this for your folder. This worked for me, the other options (commented out below) dit not.
my_options$workdir <- "C:/Users/cmagnee/Documents/cplex"
my_options$outputFolder <- "C:/Users/cmagnee/Documents/cplex"

# my_options$workdir <- "H:/COSMOS_EXIM/cplex"
# my_options$outputFolder <- "H:/COSMOS_EXIM/cplex"

# my_options$workdir <- file.path(getwd(), "cplex")
# my_options$outputFolder <- file.path(getwd(), "cplex")

my_options$solver <- "cplex" #or cbc

# max_depth <- 10
# 
# seconds_per_step <- 3600/5
# my_options$solver <- "cbc"
my_options$timelimit <- seconds_per_step
my_options$mipGAP <- 0.05

# amc computer
my_options$threads <- 19

# laptop
# my_options$threads <- 6




time_df <- data.frame()
for (patient in names(cosmos_inputs[1:4])){
  # patient <- "CLL2479_T48" # This is now automatic. Sheeeee
  sig_input <- cosmos_inputs[[patient]]$TF_scores
  metab_input <- cosmos_inputs[[patient]]$metabolomic 
  RNA_input <- cosmos_inputs[[patient]]$RNA
  
  

  # Filter inputs
  metab_input <- prepare_metab_inputs(metab_input, c("c","m"))
  # sig_input <- sig_input[abs(sig_input) > 1.5] # This is off for now
  # metab_input <- metab_input[abs(metab_input) > 1.5] # This is off for now
  metab_input <- cosmosR:::filter_input_nodes_not_in_pkn(metab_input, meta_network)
  sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)
  
  # Stp 1
  test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = meta_network[,c(1,3,2)],
                                                        signaling_data = sig_input,
                                                        metabolic_data = metab_input,
                                                        diff_expression_data = RNA_input,
                                                        maximum_network_depth = max_depth+1,
                                                        remove_unexpressed_nodes = T,
                                                        filter_tf_gene_interaction_by_optimization = T,
                                                        CARNIVAL_options = my_options)
  my_options$timelimit <- seconds_per_step
  
  #Stp 2
  test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                        CARNIVAL_options = my_options)
  
  
  formatted_res <- format_COSMOS_res(test_result_for)
  
  SIF <- formatted_res[[1]]
  ATT <- formatted_res[[2]]
  
  SIF <- SIF[which(SIF$Weight != 0),]
  
  RNA_input_df <- data.frame(Nodes = names(RNA_input), t = RNA_input)
  ATT <- merge(ATT, RNA_input_df, all.x = T)
  ATT <- ATT[ATT$AvgAct != 0,]
  
  
  write_csv(SIF, file = paste("../results/",paste(patient, "_SIF.csv",sep = ""), sep = ""))
  write_csv(ATT, file = paste("../results/",paste(patient, "_ATT.csv",sep = ""), sep = ""))
  
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
  
  #Stp 4
  my_options$timelimit <- seconds_per_step
  
  test_result_back <- run_COSMOS_metabolism_to_signaling(data = test_back,
                                                         CARNIVAL_options = my_options)
  
  formatted_res_back <- format_COSMOS_res(test_result_back)
  
  SIF_back <- formatted_res_back[[1]]
  ATT_back <- formatted_res_back[[2]]
  
  SIF_back <- SIF_back[which(SIF_back$Weight != 0),]
  
  ATT_back <- merge(ATT_back, RNA_input_df, all.x = T)
  ATT_back <- ATT_back[ATT_back$AvgAct != 0,]
  
  write_csv(SIF_back, file = paste("../results/",paste(patient, "_SIF_back.csv",sep = ""), sep = ""))
  write_csv(ATT_back, file = paste("../results/",paste(patient, "_ATT_back.csv",sep = ""), sep = ""))
  
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
  
  write_csv(SIF_full, file = paste("../results/d", paste(max_depth, paste(patient, "_SIF_full.csv", sep = ""), sep = ""), sep = ""))
  write_csv(ATT_full, file = paste("../results/d", paste(max_depth, paste(patient, "_ATT_full.csv",sep = ""), sep = ""), sep = ""))
  
  # How long did this take?
  new_row <- list(
    patient = patient,
    end_time = round(Sys.time()-start.time,2)
  )
  time_df <- rbind(time_df, new_row)
  
}
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken

