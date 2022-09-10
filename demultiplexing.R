

# Libraries  --------------------------------------------------------------
library(string)
library(tidyverse)
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

# Function to process each file: read 10x_h5, create 3 assays(RNA, protein, HTO)
createSeurat <- function(f) {
  hdf5_obj <- Read10X_h5(f) #obj contains gene expression and antibody capture data
  rn <- rownames(hdf5_obj$`Antibody Capture`)
  hashtag_names<-rn[(grep("Hashtag",rownames(hdf5_obj$`Antibody Capture`)))]
  file_name <- f 
  project_name <- substring(str_remove(file_name, gsub('.*s[0-9]', '\\1',file_name)), 32)

  #by default, the seurat object contains an assay storing RNA measurement 
  seurat_obj <- CreateSeuratObject(counts = hdf5_obj[["Gene Expression"]], project = project_name,
  min.cells = 3, min.features = 200)
  
  seurat_obj[["Protein"]] <- CreateAssayObject(head(hdf5_obj[[2]],-4)[,colnames(x=seurat_obj)])
  seurat_obj[["HTO"]] <- CreateAssayObject(tail(hdf5_obj[[2]],4)[,colnames(x=seurat_obj)])
  return(seurat_obj)
  
}
# Find all sample_filtered_feature_bc_matrix.h5 files
files <- dir( recursive=TRUE, full.names=TRUE, pattern="sample_filtered_feature_bc_matrix.h5")

# Apply the function to all files // creates multimodal seurat object 
result <- sapply(files, createSeurat)

#First perform merge, then normalize 
mSeurat <- merge(result, add.cell.ids = c("s1","s2","s3","s4", "s5","s6","s7","s8"))
#Why using both CLR and Log - normalization? 

