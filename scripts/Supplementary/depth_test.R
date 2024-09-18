################################## Libraries ################################## 
library(cosmosR)
library(readr)
library(dplyr)

################################## Importing ################################## 

T0_input_path <- "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\CosmosR\\T0_cosmos_inputs.RData"
T48_input_path <- "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\CosmosR\\T48_cosmos_inputs.RData"

data("meta_network")

meta_network <- meta_network_cleanup(meta_network)

load(T0_input_path)
load(T48_input_path)


################################## Parameters ################################## 

seconds_per_step <- 1

#In order to adapt options to users specification we can load them into a variable 
#that will then be passed to preprocess_COSMOS_signaling_to_metabolism CARNIVAL_options parameter
my_options <- default_CARNIVAL_options(solver = "cplex")

#Here the user should provide a path to its CPLEX executable

# my_options$solverPath <- "~/Documents/cplex" #or cbc solver executable

my_options$solverPath <- "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64"

# Change this for your folder. This worked for me, the other options (commented out below) did not.
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


# Initialize an empty dataframe to store the results
filter_df <- data.frame(patient = character(),
                        n_sig_input = integer(),
                        f_sig = numeric(),
                        n_metab_input = integer(),
                        f_metab = numeric(),
                        stringsAsFactors = FALSE)

filter_s <- 2.0  # Start filter_s 
filter_m <- 0.7  # Start filter_m 
cosmos_inputs <- T0_cosmos_inputs # T0 or T48
for (i in 1:2){
  cosmos_inputs <- list(T0_cosmos_inputs, T48_cosmos_inputs)[[i]]
  for (patient in names(cosmos_inputs)) {
    max_depth <- 2  # Initialize max_depth for each patient
    while (max_depth <= 10) {  # Continue while max_depth <= 10
      skip_to_next <- FALSE
      
      tryCatch({
        # Extract patient data
        sig_input <- cosmos_inputs[[patient]]$TF_scores
        metab_input <- cosmos_inputs[[patient]]$metabolomic
        RNA_input <- cosmos_inputs[[patient]]$RNA
        
        # Prepare metabolomic inputs
        metab_input <- prepare_metab_inputs(metab_input, c("c", "m"))
        
        # Filter inputs based on thresholds
        sig_input <- sig_input[abs(sig_input) > filter_s]
        metab_input <- metab_input[abs(metab_input) > filter_m]
        
        # Filter nodes not present in the meta_network
        metab_input <- cosmosR:::filter_input_nodes_not_in_pkn(metab_input, meta_network)
        sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)
        
        # Preprocess signaling to metabolism
        test_for <- preprocess_COSMOS_signaling_to_metabolism(
          meta_network = meta_network[, c(1, 3, 2)],
          signaling_data = sig_input,
          metabolic_data = metab_input,
          diff_expression_data = RNA_input,
          maximum_network_depth = max_depth,
          remove_unexpressed_nodes = TRUE,
          filter_tf_gene_interaction_by_optimization = TRUE,
          CARNIVAL_options = my_options
        )
        
        # Preprocess metabolism to signaling
        test_back <- preprocess_COSMOS_metabolism_to_signaling(
          meta_network = meta_network[, c(1, 3, 2)],
          signaling_data = sig_input,
          metabolic_data = metab_input,
          diff_expression_data = RNA_input,
          maximum_network_depth = max_depth,
          remove_unexpressed_nodes = TRUE,
          filter_tf_gene_interaction_by_optimization = TRUE,
          CARNIVAL_options = my_options
        )
        
        # If no error occurs, add the current row to filter_df
        new_row <- data.frame(
          patient = patient,
          max_depth = max_depth,
          n_sig_input = length(sig_input),
          n_metab_input = length(metab_input),
          stringsAsFactors = FALSE
        )
        
        filter_df <- rbind(filter_df, new_row)
        print(filter_df)
        
        skip_to_next <- TRUE  # Indicate successful processing
        
      }, error = function(e) {
        # Error handling: move to the next depth level
        message(paste("Error in processing patient:", patient, "-", e$message))
        max_depth <<- max_depth + 1  # Increment max_depth
        
      })
      
      # Break the while loop and move to the next patient if successful
      if (skip_to_next) {
        break  # Exit the while loop for the current patient
      }
    }
  }
}

