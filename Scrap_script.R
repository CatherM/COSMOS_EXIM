library(readxl)
library(readr)
library(biomaRt)
library(dplyr)
library(reshape2)
library(pheatmap)
library(vsn)

library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)


## Gene expression
rawcounts_path <- "L:\\basic\\divg\\EXIM\\ImmunoHematology\\Cathy Magnée\\Data\\CD3_HDvsCLL\\Transcriptomic\\EC024_raw counts.csv"
# rsem_count.path <- "L:/basic/divg/EXIM/ImmunoHematology/Elena Camerini/Bioinformatics/EC024_CJ_RNAseq/Raw data from alignment/rsem_count"

# Loading data, changing first column into rownames and removing this column.
RNA_rawcounts <- as.data.frame(read_csv(rawcounts_path))

# Make rownames, then remove them (only to put them back later)
rownames(RNA_rawcounts) <- RNA_rawcounts$...1
RNA_rawcounts <- RNA_rawcounts[,c(-1)]


# Order them so Day 0 and day 2 are separate. 
RNA <- RNA_rawcounts[,c(7,9,6,2,5,14,15,16,17,18,3,8,4,1,10,11,12,13)] 
RNA$GENE_ID <- rownames(RNA)

# Remove redundant variables
rm(rawcounts_path)

day0_RNA_raw <- RNA[,c(1:10)]
day2_RNA_raw <- RNA[,c(11:18)]


myCPM <- cpm(day0_RNA_raw)
tresh <- myCPM > 0.5

table(rowSums(tresh))
keep <- rowSums(tresh) >= 4
counts.keep <- day0_RNA_raw[keep,]
summary(keep)
dim(counts.keep)

