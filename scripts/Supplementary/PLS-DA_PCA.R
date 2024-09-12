# PLS-DA_PCA.R
# Script goal:
# to find out if the metabolic data still contains data that is able to differentiate between CLL and HD, after preprocessing steps

# Libraries:
library(readr)
library(ropls) #install via BiocManager
library(mixOmics)
library(pheatmap)

# Data:
path_T0_metabolomic_data <- "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\Metabolomic\\T0_metabs_hmdb_scaled.csv"
T0_metabolomic <- as.data.frame(read_csv(path_T0_metabolomic_data, show_col_types = FALSE))
names(T0_metabolomic)[names(T0_metabolomic) == "...1"] <- "metabolite"
rownames(T0_metabolomic) <- T0_metabolomic$metabolite
T0_metabolomic <- subset(T0_metabolomic, select=-c(metabolite))
T0_metabolomic_wide <- data.frame(t(T0_metabolomic))

path_T48_metabolomic_data <- "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\Metabolomic\\T48_metabs_hmdb_scaled.csv"
T48_metabolomic <- as.data.frame(read_csv(path_T48_metabolomic_data, show_col_types = FALSE))
names(T48_metabolomic)[names(T48_metabolomic) == "...1"] <- "metabolite"
rownames(T48_metabolomic) <- T48_metabolomic$metabolite
T48_metabolomic <- subset(T48_metabolomic, select=-c(metabolite))
T48_metabolomic_wide <- data.frame(t(T48_metabolomic))

rm(path_T0_metabolomic_data) 
rm(path_T48_metabolomic_data) 

#########################
#########################
#########################
test_df <- T48_metabolomic_wide
metadata_df <- data.frame()

# Create metadata df
for (patient in rownames(test_df)){
  if (startsWith(patient, "HD")){
    new_row <- list(
      Sample = patient, 
      Disease_state = "HD")
    metadata_df <- rbind(metadata_df, new_row)
  } else {
    new_row <- list(
      Sample = patient, 
      Disease_state = "CLL")
    metadata_df <- rbind(metadata_df, new_row)
  }
}
rm(patient)
rm(new_row)

rownames(metadata_df) <- metadata_df$Sample
metadata_df$Disease_state <-  as.factor(metadata_df$Disease_state)
Disease_state_fc <- metadata_df[, "Disease_state"]


# color_mapping <- c("CLL" = "red", "HD" = "blue")

# Assign colors based on the Disease_state_fc factor
# sample_colors <- color_mapping[metadata_df$Disease_state]

# Create pca graphs
test_df.pca <- opls(test_df)

# Plot pca graphs with custom colors for disease state
plot(test_df.pca,
     typeVc = "x-score",
     parAsColFcVn = Disease_state_fc,
     parPaletteVc = c("red", "blue"))

# plsda graph
# sacurine.plsda <- opls(test_df, Disease_state_fc, parPaletteVc = c("red", "blue"))


# Heatmmap based on VIP score
# Install and load necessary packages
# Create a vector indicating the group labels
# group <- factor(c(rep("CLL", 4), rep("HD", 5)))

# Perform PLS-DA
# plsda_res <- plsda(T0_metabolomic_wide, group, ncomp = 2)

# Extract VIP scores for each metabolite
# vip_scores <- vip(plsda_res)

# Sort VIP scores in decreasing order and select top 50
# top_50_indices <- order(vip_scores, decreasing = TRUE)[1:50] ########## MISTAKE

# Subset the original data to include only these top 50 metabolites
# T0_metabolomic_top50 <- T0_metabolomic[top_50_indices, ]
# 
# pheatmap(as.matrix(complete.cases(T0_metabolomic_top50)), 
#          scale = "row", 
#          clustering_distance_rows = "euclidean", 
#          clustering_distance_cols = "euclidean", 
#          clustering_method = "complete", 
#          annotation_col = annotation_col)
# 
# pheatmap(as.matrix(complete.cases(T0_metabolomic_top50)), 
#          scale = "row",  # Scale rows (metabolites) to have mean=0 and variance=1
#          clustering_distance_rows = "euclidean", 
#          clustering_distance_cols = "euclidean", 
#          clustering_method = "complete", 
#          annotation_col = data.frame(Group = group))
