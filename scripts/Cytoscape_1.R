# Cytoscape.R
# First Scrap Cytoscape Script. 
# Script goal: Exploring possibilities of automating Cytoscape via RStudio.

# NOTE: Make sure Cytoscape is open on your PC before going

# Library start:
library(RCy3)
library(readr)
library(stringr)
# library(colorRamps)
# library(STRINGdb)

# OPEN CYTOSCAPE BEFORE CONTINUING
# Check Connection to Cytoscape:
cytoscapePing()
cytoscapeVersionInfo()
style.name = "empty"
# If the above gives errors, please restart R and cytoscape and try again. 

# Import dataset so you can loop through the patient names:
# Import datasets:

# T48_input_path <- "ourdata/cosmos/T48_cosmos_inputs.RData" # Laptop only
# load(T48_input_path)
# T48_cosmos_inputs <- day2_cosmos_inputs # Laptop only
# rm(day2_cosmos_inputs)
# cosmos_inputs <- T48_cosmos_inputs


# COSMOS variables (What network do you want to see / make?). Folder path
max_depth <- 3
filter <- "s15m05"
stim_point <- "T48"

# AMC PATH BLOCK
core_path <- "L:/basic/divg/EXIM/ImmunoHematology/Cathy Magnée/Data/CD3_HDvsCLL/CosmosR"
result_path <- paste("/Networks/d", max_depth, sep = "") # Add depth folder
result_path <- paste(result_path, paste("/f", filter, sep = ""), sep = "") # Add filter folder
result_path <- paste(result_path, paste("/", stim_point, sep = ""), sep = "") # Add stim/unstim point folder
result_path <- paste(result_path, "/csvs/", sep = "") # Add csv folder
data_path <- paste(core_path, result_path, sep = "")

# Laptop path block:

# core_path <- "C:/Users/Cathy.LAPTOP-SDFOSKVI/Documents/School_files/SysBio_RP/COSMOS_EXIM/results/processed_results"
# result_path <- paste("/d", max_depth, sep = "") # Add depth folder
# data_path <- paste(core_path, result_path, sep = "")


# Create Merged_Network ATT dfs
Merged_ATT_CLL <- data.frame()
Merged_ATT_HD <- data.frame()

# Create Merged_Network SIF dfs
Merged_SIF_CLL <- data.frame()
Merged_SIF_HD <- data.frame()

all_files <- list.files(path = data_path, pattern = "*_SIF_processed.csv", full.names = TRUE)
for (files in all_files){
  patient <- str_extract(files, "[^/]+(?=_SIF)")
  # Import processed SIF and ATT files for patient
  SIF_path <- paste(data_path, paste(patient, "_SIF_processed.csv", sep = ""), sep = "")
  SIF <- read_csv(SIF_path, show_col_types = F)
  SIF <- unique(SIF)
  ATT_path <- paste(data_path, paste(patient, "_ATT_processed.csv",sep = ""), sep = "")
  ATT <- read_csv(ATT_path, show_col_types = F)
  ATT <- unique(ATT)
  
  if(substr(patient, 0, 3)== "CLL"){ 
    for (i in 1:nrow(SIF)){
      Edge <- SIF[i,]
      Merged_SIF_CLL <- rbind(Merged_SIF_CLL, Edge)
    }
    # Only for CLL patients
    for (Noderow in unique(ATT$Nodes)){
      # For every node in the ATT 
      if (Noderow %in% Merged_ATT_CLL$Nodes){
        Merged_ATT_CLL$Measured[Merged_ATT_CLL$Nodes == Noderow] <- 
          Merged_ATT_CLL$Measured[Merged_ATT_CLL$Nodes == Noderow] + 1
        Merged_ATT_CLL$Activity[Merged_ATT_CLL$Nodes == Noderow] <- 
          Merged_ATT_CLL$Activity[Merged_ATT_CLL$Nodes == Noderow] + mean(ATT$Activity[ATT$Nodes == Noderow])
      } else {
        new_row <- list(
          Nodes = Noderow,
          Measured = 1,
          measured = mean(ATT$measured[ATT$Nodes == Noderow]),
          MoleculeType = ATT$MoleculeType[ATT$Nodes == Noderow][1],
          Activity = mean(ATT$Activity[ATT$Nodes == Noderow])
        )
        Merged_ATT_CLL <- rbind(Merged_ATT_CLL, new_row)
      }
    }

    
  } else {
    for(Noderow in unique(ATT$Nodes)){
      # Only for HD patients
      if (Noderow %in% Merged_ATT_HD$Nodes){
        Merged_ATT_HD$Measured[Merged_ATT_HD$Nodes == Noderow] <- 
          Merged_ATT_HD$Measured[Merged_ATT_HD$Nodes == Noderow] + 1
        Merged_ATT_HD$Activity[Merged_ATT_HD$Nodes == Noderow] <- 
          Merged_ATT_HD$Activity[Merged_ATT_HD$Nodes == Noderow] + mean(ATT$Activity[ATT$Nodes == Noderow])
      } else {
        new_row <- list(
          Nodes = Noderow,
          Measured = 1,
          measured = mean(ATT$measured[ATT$Nodes == Noderow]),
          MoleculeType = ATT$MoleculeType[ATT$Nodes == Noderow][1],
          Activity = mean(ATT$Activity[ATT$Nodes == Noderow])
        )
        Merged_ATT_HD <- rbind(Merged_ATT_HD, new_row)
      }
    }
    for (i in 1:nrow(SIF)){
      Edge <- SIF[i,]
      Merged_SIF_HD <- rbind(Merged_SIF_HD, Edge)
    }
  }
}

# Only keep unique edges
Merged_SIF_HD <- unique(Merged_SIF_HD) 
Merged_SIF_CLL <- unique(Merged_SIF_CLL)

# Only keep Nodes that are present in more than 1 network
Merged_ATT_HD <- filter(Merged_ATT_HD, Measured > 1)
Merged_ATT_CLL <- filter(Merged_ATT_CLL, Measured > 1)

# Only keep edges for which both nodes are present in the network
Merged_SIF_CLL <- filter(Merged_SIF_CLL, Node1 %in% Merged_ATT_CLL$Nodes)
Merged_SIF_CLL <- filter(Merged_SIF_CLL, Node2 %in% Merged_ATT_CLL$Nodes)

Merged_SIF_HD <- filter(Merged_SIF_HD, Node1 %in% Merged_ATT_HD$Nodes)
Merged_SIF_HD <- filter(Merged_SIF_HD, Node2 %in% Merged_ATT_HD$Nodes)


###################################################################################### 
############################ Create network function #################################
###################################################################################### 

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
  ATT_string_interaction_cmd <- paste('string protein query taxonID=9606 cutoff=0.9 query=', paste(nodes$id, collapse=","),'"',sep="")
  commandsGET(ATT_string_interaction_cmd )
  string_nodes <- getTableColumns(table = "node", columns = c("display name", "target::family"), network = "STRING network")
  colnames(string_nodes) <- c("id", "moleculeType")
  string_nodes <- merge(nodes, string_nodes, by = "id", all.x = T)
  nodes$moleculeType <- ifelse(is.na(string_nodes$moleculeType.x),
                                     string_nodes$moleculeType.y,
                                     string_nodes$moleculeType.x)
  nodes$moleculeType[is.na(nodes$moleculeType)] <- "Other"
  
  
  closeSession(save.before.closing = F)

  
  createNetworkFromDataFrames(nodes, edges, title = MergedGroup, collection = MergedGroup) # if(substr(MergedGroup, 0, 3)== "CLL"){"CLL"} else { "HD"})
  # You should see a network rn frfr
  
  # Setting style for network OR creating it first, then setting the style. If else so we do not remake this style over and over.
  if (style.name == "COSMOS_Style"){
    setVisualStyle(style.name)
  } else {
    style.name = "COSMOS_Style"
    defaults <- list(NODE_SHAPE="ROUND_RECTANGLE",
                     NODE_SIZE=50,
                     EDGE_TRANSPARENCY=255,
                     NODE_LABEL_POSITION="S,N,c,0.00,0.00")
    nodeLabels <- mapVisualProperty('node label','id','p')
    nodeFills <- mapVisualProperty('node fill color',
                                   'Activity',
                                   'c',
                                   c(-1,1),
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
  saveSession(filename = paste(data_path, MergedGroup, sep = "Merged_Network_"))
  
}


# Creating networks here, using Create_Network function: 
Create_Network(ATT = Merged_ATT_CLL, SIF= Merged_SIF_CLL, MergedGroup = "CLL")

Create_Network(ATT = Merged_ATT_HD, SIF= Merged_SIF_HD, MergedGroup = "HD")


####################### Any individual Network here: #############################

patient = "CLL2105_T48" # Fill it in
SIF_path <- paste(data_path, paste(patient, "_SIF_processed.csv", sep = ""), sep = "")
SIF <- read_csv(SIF_path, show_col_types = F)
ATT_path <- paste(data_path, paste(patient, "_ATT_processed.csv",sep = ""), sep = "")
ATT <- read_csv(ATT_path, show_col_types = F)

Create_Network(ATT = ATT, SIF= SIF, MergedGroup = patient)
