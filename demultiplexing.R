

# Libraries  --------------------------------------------------------------
library(stringr)
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
vector_result <- unlist(result)
mSeurat <- merge(vector_result[[1]],y = vector_result[2:8], add.cell.ids = c("s1","s2","s3","s4", "s5","s6","s7","s8"))
#By default, merge() will combine the Seurat objects based on the raw count matrices,
#erasing any previously normalized and scaled data matrices -> normalizing before merging does not make sense 

#Why using both CLR and Log - normalization? So in the end, which one? 

# QC  ---------------------------------------------------------------------
# % MT reads
mSeurat[["percent.mt"]] <- PercentageFeatureSet(mSeurat, pattern = "^MT-")
head(mSeurat@meta.data)
View(mSeurat@meta.data)

VlnPlot(mSeurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#looking at features together 
FeatureScatter(mSeurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')


# Filtering ---------------------------------------------------------------
mSeurat <- subset(mSeurat, subset = percent.mt < 10 & nFeature_RNA > 500)
dim(mSeurat)


# Normalization -----------------------------------------------------------
mSeurat <- NormalizeData(mSeurat)

# Changing the format of cell names for later
cellnames <- paste0(sapply(strsplit(colnames(mSeurat), "_|-"),"[[",2),
                    "-",
                    sapply(strsplit(colnames(mSeurat), "_|-"),"[[",1))
print(head(colnames(mSeurat)))
print(head(cellnames))
# workaround
cellnames <- gsub("s", "", cellnames)

mSeurat <- RenameCells(mSeurat, new.names = cellnames)


# Finding Variable Features -----------------------------------------------
mSeurat <- FindVariableFeatures(mSeurat, selection.method = "vst", nfeatures = 5000)
top10 <- head(VariableFeatures(mSeurat), 10)
plot1 <- VariableFeaturePlot(mSeurat)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Hashtag assign -----------------------------------------------------------------

mSeurat@meta.data$hash.ID<-apply(
  as.matrix(mSeurat@assays$HTO@counts), 2, 
  function (col) rownames(as.matrix(mSeurat@assays$HTO@counts))[which.max(col)])
# if most highly hashtag does not make up at least 75%, classify as doublet
proportion <- apply(as.matrix(mSeurat@assays$HTO@counts), 2, max) / 
  colSums(as.matrix(mSeurat@assays$HTO@counts))
# View(as.matrix(mSeurat@assays$HTO@counts)[,which(proportion < 0.8)])
mSeurat$hash.ID[which(proportion < 0.8)] <- "doublets"
mSeurat <- subset(mSeurat, hash.ID == "doublets", invert = T) # and remove
remove(proportion)

#mSeurat <- HTODemux(mSeurat, kfunc = "kmeans", init = 3, positive.quantile = 0.95)
#table(mSeurat@meta.data$hash.ID)
#HTOHeatmap(mSeurat)
#table(mDemux$HTO_classification.global)


 # Group cells based on the max HTO signal
 Idents(mSeurat) <- "HTO_maxID"
 RidgePlot(mSeurat, assay = "HTO", features = rownames(mSeurat[["HTO"]])[1:4], ncol = 2)
 
 mSeurat$sample <- paste(mSeurat$orig.ident, mSeurat$hash.ID)
 Idents(mSeurat) <- "sample"

 mSeurat <- subset(
   mSeurat, idents = c("s1 Hashtag1","s1 Hashtag2",
                       "s2 Hashtag1","s2 Hashtag2",
                       "s3 Hashtag3","s3 Hashtag4",
                       "s4 Hashtag3","s4 Hashtag4",
                       "s5 Hashtag1","s5 Hashtag2",
                       "s6 Hashtag1","s6 Hashtag2",
                       "s7 Hashtag3","s7 Hashtag4",
                       "s8 Hashtag3","s8 Hashtag4"))
 print(paste("Number of cells:", ncol(mSeurat)))

 mSeurat$sample <- stringi::stri_replace_all_regex(
   mSeurat$sample,
   sort(unique(mSeurat$sample)),
   c("Dnr58_Tonsil_CXCR5positive", "Dnr58_Blood_CXCR5positive", 
     "Dnr58_Tonsil_CXCR5negative", "Dnr58_Blood_CXCR5negative", 
     
     "Dnr69_Tonsil_CXCR5positive", "Dnr69_Blood_CXCR5positive", 
     "Dnr69_Tonsil_CXCR5negative", "Dnr69_Blood_CXCR5negative", 
     
     "Dnr72_Tonsil_CXCR5positive", "Dnr72_Blood_CXCR5positive", 
     "Dnr72_Tonsil_CXCR5negative", "Dnr72_Blood_CXCR5negative", 
     
     "Dnr83_Tonsil_CXCR5positive", "Dnr83_Blood_CXCR5positive", 
     "Dnr83_Tonsil_CXCR5negative", "Dnr83_Blood_CXCR5negative"),
   vectorize=FALSE)
 Idents(mSeurat) <- "sample"
 table(mSeurat$sample)
 