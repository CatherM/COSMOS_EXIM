# Cytoscape.R
# First Scrap Cytoscape Script. 
# Script goal: Exploring possibilities of automating Cytoscape via RStudio.

# NOTE: Make sure Cytoscape is open on your PC before going

# Library start:
library(RCy3)
library(readr)

# Path variables (What network do you want to see / make?)
max_depth <- 4
filter <- 125
stim_point <- "T0"



core_path <- "L:/basic/divg/EXIM/ImmunoHematology/Cathy Magnée/Data/CD3_HDvsCLL/CosmosR"

result_path <- paste("/Networks/d", max_depth, sep = "") # Add depth folder
result_path <- paste(result_path, paste("/f", filter, sep = ""), sep = "") # Add filter folder
result_path <- paste(result_path, paste("/", stim_point, sep = ""), sep = "") # Add stim/unstim point folder
result_path <- paste(result_path, "/csvs/", sep = "") # Add csv folder

data_path <- paste(core_path, result_path, sep = "")

# For now, we will do it per sample. But later we might try to loop through the datasets T0 and T48

patient <- "HD91_T0"

SIF_path <- paste(data_path, paste(patient, "_SIF_full.csv", sep = ""), sep = "")
SIF <- read_csv(SIF_path)

ATT_path <- paste(data_path, paste(patient, "_ATT_full.csv",sep = ""), sep = "")
ATT <- read_csv(ATT_path)

# data <- SIF

nodes <- data.frame(id = ATT$Nodes,
                    Activity = ATT$`Activity_ HD91_T0`,
                    measured = ATT$measured,
                    stringsAsFactors = F)
edges <- data.frame(source = SIF$Node1,
                    target = SIF$Node2,
                    interaction = SIF$Sign,
                    weight = SIF$Weight,
                    stringsAsFactors = F)

createNetworkFromDataFrames(nodes, edges, title = "HD91", collection = "T0")

