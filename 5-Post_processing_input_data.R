## Dit is een kladblad voor het post-preprocessen van het netwerk. 
## Goal: Vinden welke nodes de inputs zijn

#### This is how you find out which of the Nodes are input

# RNA
intersect(names(cosmos_inputs$`RPMI-8226`$RNA), SIF_full$Node1)
intersect(names(cosmos_inputs$`RPMI-8226`$RNA), SIF_full$Node2)

# TFs
intersect(names(cosmos_inputs$`RPMI-8226`$TF_scores), SIF_full$Node1)
intersect(names(cosmos_inputs$`RPMI-8226`$TF_scores), SIF_full$Node2)

# Metabolic does not work as they are HMDB 
# Still need something for this
# intersect(names(cosmos_inputs$`RPMI-8226`$metabolomic), SIF_full$Node1)


#### Finding if they are up- or downregulated
# Change this to any of the above 4-6 


xy <- intersect(names(cosmos_inputs$`RPMI-8226`$TF_scores), SIF_full$Node2)


# Dont forget to change data type here (TF_scores or RNA)
for (node in xy){
  print(paste0(node,": " ,cosmos_inputs$`RPMI-8226`$TF_scores[[node]]))
}


### Extra try to change metabs to their non-hmdb names
metabolite_matching <- as.data.frame(read_csv("Factor_COSMOS/support/metabolite_matching.csv"))


# Function to process nodes and return corresponding HMDB IDs
get_hmdb_ids <- function(SIF_full, metabolite_matching) {
  # Combine Node1 and Node2 into one vector
  all_nodes <- unique(c(SIF_full$Node1, SIF_full$Node2))
  
  # Filter nodes containing "Metab_"
  metab_nodes <- all_nodes[grep("Metab_", all_nodes)]
  
  # Remove "Metab_" and everything after the last "_"
  processed_nodes <- sub("Metab_", "", metab_nodes)
  processed_nodes <- sub("_[^_]*$", "", processed_nodes)
  
  # Approximate matching of processed nodes with metabolites
  matched_indices <- amatch(processed_nodes, metabolite_matching$metabolite, maxDist = 10)
  
  # Retrieve corresponding HMDB IDs
  matched_hmdb <- metabolite_matching$hmdb[matched_indices]
  
  return(matched_hmdb)
}

hmdb_ids <- get_hmdb_ids(SIF_full, metabolite_matching)
print(hmdb_ids)
# subset(SIF_full$Node1, grepl("Metabs_"))
# for (node in SIF_full$Node1){
#   if node contains()
# }