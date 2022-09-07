

# Libraries  --------------------------------------------------------------
library(string)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(patchwork)
library(gtools)
library(stringr)
library(djvdj)
library(ggrepel)
library(ggpubr)
library(tidyr)
library(EnhancedVolcano)
library(RColorBrewer)
library(msigdbr)
library(fgsea)
library(gridExtra)



# Data import -----------------------------------------------------------



# Create Seurat Objects with raw (non - normalized data)  -----------------

# Function to process each file
createSeurat <- function(f) {
  hdf5_obj <- Read10X_h5(f) #obj contains gene expression and antibody capture data 
  project_name <- substring(str_remove(files, gsub('.*s[0-9]', '\\1',f)), 6)
  seurat_obj <- CreateSeuratObject(counts = hdf5_obj, project = project_name,
  min.cells = 3, min.features = 200)
}

# Find all sample_filtered_feature_bc_matrix.h5 files
files <- dir( recursive=TRUE, full.names=TRUE, pattern="sample_filtered_feature_bc_matrix.h5")

# Apply the function to all files
result <- sapply(files, createSeurat)
