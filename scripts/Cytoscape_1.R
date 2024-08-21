# Cytoscape.R
# First Scrap Cytoscape Script. 
# Script goal: Exploring possibilities of automating Cytoscape via RStudio.

# NOTE: Make sure Cytoscape is open on your PC before going

# Library start:
library(RCy3)
library(readr)
library(colorRamps)
library(STRINGdb)

# OPEN CYTOSCAPE BEFORE CONTINUING
# Check Connection to Cytoscape:
cytoscapePing()
cytoscapeVersionInfo()
style.name = "empty"
# If the above gives errors, please restart R and cytoscape and try again. 

# Import dataset so you can loop through the patient names:
# Import datasets:
T48_input_path <- "ourdata/cosmos/T48_cosmos_inputs.RData" # Laptop only
load(T48_input_path)
T48_cosmos_inputs <- day2_cosmos_inputs # Laptop only
rm(day2_cosmos_inputs)
cosmos_inputs <- T48_cosmos_inputs


# COSMOS variables (What network do you want to see / make?). Folder path
max_depth <- 2
filter <- 125
stim_point <- "T0"

# AMC PATH BLOCK
core_path <- "L:/basic/divg/EXIM/ImmunoHematology/Cathy Magnée/Data/CD3_HDvsCLL/CosmosR"
result_path <- paste("/Networks/d", max_depth, sep = "") # Add depth folder
result_path <- paste(result_path, paste("/f", filter, sep = ""), sep = "") # Add filter folder
result_path <- paste(result_path, paste("/", stim_point, sep = ""), sep = "") # Add stim/unstim point folder
result_path <- paste(result_path, "/csvs/", sep = "") # Add csv folder
data_path <- paste(core_path, result_path, sep = "")

# Laptop path block:

core_path <- "C:/Users/Cathy.LAPTOP-SDFOSKVI/Documents/School_files/SysBio_RP/COSMOS_EXIM/results/processed_results"
result_path <- paste("/d", max_depth, sep = "") # Add depth folder
data_path <- paste(core_path, result_path, sep = "")


# Create Merged_Network dfs
Merged_Network_CLL <- data.frame()
Merged_Network_HD <- data.frame()
patient <- "CLL2648_T48"
for (patient in names(cosmos_inputs)){
  
  # Import processed SIF and ATT files for patient
  SIF_path <- paste(data_path, paste(patient, "_SIF_processed.csv", sep = ""), sep = "")
  SIF <- read_csv(SIF_path)
  SIF <- unique(SIF)
  ATT_path <- paste(data_path, paste(patient, "_ATT_processed.csv",sep = ""), sep = "")
  ATT <- read_csv(ATT_path)
  ATT <- unique(ATT)
  
  
  if(substr(patient, 0, 3)== "CLL"){
    for (Noderow in unique(ATT$Nodes)){
      if (Noderow in Merged_Network_CLL$Nodes){
        continue
      } else {
        new_row <- list(
          Nodes = Noderow
          Activity = ATT$Activity[Noderow]
        )
        Merged_Network_CLL <- rbind(Merged_Network_CLL)
      }
    }
    
  } else { 
      "HD"
    }
  
  # Translate this to Node Table and Edges Table for in Cytoscape
  nodes <- data.frame(id = ATT$Nodes,
                      Activity = ATT$`Activity`,
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
  
  createNetworkFromDataFrames(nodes, edges, title = patient, collection = if(substr(patient, 0, 3)== "CLL"){"CLL"} else { "HD"})
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
               'GCPR', 'TF-Epigenetic', 
               'Nuclear Receptor', 'Transporter')
  shapes <- c ('ELLIPSE', 'TRIANGLE', 
               'DIAMOND', 'PARALLELOGRAM', 
               'RECTANGLE', 'VEE', 
               'HEXAGON', 'DIAMOND',
               'OCTAGON', 'ROUND RECTANGLE')
  setNodeShapeMapping (column, values, shapes, style.name = style.name)
  
  getLayoutNames()
  getLayoutPropertyNames() # look up what these mean and fill in list below.
  setLayoutProperties("degree-circle", properties.list = list())
}