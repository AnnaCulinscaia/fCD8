

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
library(ggplot2)




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
 
 # Add patient and timepoint by splitting from sample metadata
 mSeurat$Patient <- t(as.data.frame(strsplit(mSeurat$sample, "_")))[,1]
 mSeurat$Tissue <- t(as.data.frame(strsplit(mSeurat$sample, "_")))[,2]
 mSeurat$Group <- t(as.data.frame(strsplit(mSeurat$sample, "_")))[,3]

 locations <- c("~/22-08-ORBCXCR5/cellranger/e19s1_multi-part/vdj_t/",
                "~/22-08-ORBCXCR5/cellranger/e19s2_multi-part/vdj_t/",
                "~/22-08-ORBCXCR5/cellranger/e19s3_multi-part/vdj_t/",
                "~/22-08-ORBCXCR5/cellranger/e19s4_multi-part/vdj_t/",
                "~/22-08-ORBCXCR5/cellranger/e19s5_multi-part/vdj_t/",
                "~/22-08-ORBCXCR5/cellranger/e19s6_multi-part/vdj_t/",
                "~/22-08-ORBCXCR5/cellranger/e19s7_multi-part/vdj_t/",
                "~/22-08-ORBCXCR5/cellranger/e19s8_multi-part/vdj_t/")
 suffices <- c("1", "2", "3", "4", 
               "5", "6", "7", "8")
 
 output <- list()
 original_meta <- mSeurat@meta.data
 for (i in 1:length(locations)){
   # get the cells belonging to the current set
   positions <- grep(suffices[i], colnames(mSeurat))
   # temporarily adjust to make cell name suffix the same
   TCR <- mSeurat[,positions]
   cellnames <- colnames(TCR)
   TCR <- RenameCells(TCR, new.names = gsub(suffices[i], "1", cellnames))
   # import the tcr
   TCR <- import_vdj(TCR, vdj_dir = locations[i])
   # change the cell name back and store the data
   output[[i]] <- RenameCells(TCR, new.names = cellnames)@meta.data
 }
 
 meta_w_TCR <- do.call(rbind, output)
 # make the order of the cells the same
 meta_w_TCR <- meta_w_TCR[rownames(mSeurat@meta.data),]
 mSeurat@meta.data <- meta_w_TCR # add the collated information
 mSeurat@meta.data$clonotype_id <- NULL # No longer useful because id is repeated between sets
 remove(meta_w_TCR, i, suffices, locations, cellnames, TCR, output, original_meta, positions)
 
 ###########################################
 # For cases with three or 4 chains, only take the most highly expressed TRA TRB
 
 # Perform clean up and filtering of TCR metadata information.
 # Take only the most highly expressed TCR in instances of >2 chains
 cleanMeta <- mSeurat@meta.data
 col_to_edit <- c("v_gene","d_gene","j_gene", "c_gene", "chains","cdr3",
                  "cdr3_nt", "reads", "umis", "productive", "full_length")
 
 # Looks like there is not the case of chains == TRA;TRA or TRB;TRB
 # Which makes things easier for me
 rows_to_edit <- which(cleanMeta$n_chains > 2)
 for (r in rows_to_edit){
   expressions <- as.numeric(paste0(
     unlist(strsplit(cleanMeta[r, "umis"], ";")),
     ".",
     unlist(strsplit(cleanMeta[r, "reads"], ";"))))
   
   is.TRA <- which(unlist(strsplit(cleanMeta[r, "chains"], ";")) == "TRA")
   is.TRB <- which(unlist(strsplit(cleanMeta[r, "chains"], ";")) == "TRB")
   to.keep <- c(which(expressions == max(expressions[is.TRA])),
                which(expressions == max(expressions[is.TRB])))
   
   for (x in col_to_edit){
     cleanMeta[r,x] <- paste(unlist(strsplit(cleanMeta[r,x], ";"))[to.keep],
                             collapse = ";")
   }
 }
 mSeurat@meta.data <- cleanMeta
 remove(cleanMeta, rows_to_edit, col_to_edit, r, expressions,
        is.TRA, is.TRB, to.keep, x)

 