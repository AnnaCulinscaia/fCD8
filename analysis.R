
# make clustering and umap
mSeurat <- ScaleData(mSeurat, features = rownames(mSeurat))
mSeurat <- RunPCA(mSeurat)
ElbowPlot(mSeurat)
mSeurat <- FindNeighbors(mSeurat, dims = 1:15)
mSeurat <- FindClusters(mSeurat, resolution = 0.2)
mSeurat <- RunUMAP(mSeurat, dims = 1:15, umap.method = "uwot")
saveRDS(mSeurat, paste0("./objects/mSeurat_umap",
                       format(Sys.Date(), format="%y%m%d"),
                       ".rds"))


DimPlot(mSeurat) +
  DimPlot(mSeurat, group.by = "Tissue", shuffle = T) +
  DimPlot(mSeurat, group.by = "Group", shuffle = T) +
  DimPlot(mSeurat, group.by = "Patient", shuffle = T)



Idents(mSeurat) <- "seurat_clusters"
markers <- FindAllMarkers(mSeurat,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.75, # higher to save time
                          densify = T)
markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(mSeurat, features = c(top10$gene)) + NoLegend()
DotPlot(mSeurat, features = rev(top10$gene[!duplicated(top10$gene)]),
        group.by = "seurat_clusters") + RotatedAxis()

# add cluster names
clusters_named <- c("Cytotoxic","SomeCD69+CD103+1","SomeCD69+CD103+2",
                    "Some CD127+CD62L+","CD3 negative","CD95-FOS(high)JUNB(high) (some)",
                    "KLRB1+","CD69+CD103+CD71+ (protein)","RORA(high)RUNX2(high)")
names(clusters_named) <- levels(mSeurat)
mSeurat <- RenameIdents(mSeurat, clusters_named)
mSeurat[["clusters_named"]] <- Idents(object = mSeurat)


DimPlot(mSeurat, label = T) + NoLegend()

# cluster 4 is the non T cells
FeaturePlot(mSeurat, c("CD3D", "CD3E", "CD3G", "percent.mt"))
FeaturePlot(mSeurat, c("CD8A", "CD8B", "CD4", "CD5"))
FeaturePlot(mSeurat, "nCount_RNA",max.cutoff = c(7500))
FeaturePlot(mSeurat, "nFeature_RNA", max.cutoff = c(3000))
FeaturePlot(mSeurat, c("PTPRC", "CD14","LYZ","MS4A1"))
FeaturePlot(mSeurat, c("FCGR3A", "CST3", "CD68","TRAC"))
FeaturePlot(mSeurat, c("MZB1", "MKI67", "KIT"))


Idents(mSeurat) <- "seurat_clusters"
nonT <- FindMarkers(mSeurat, ident.1 = 4, logfc.threshold = 1)

# GSEA
#Retrieve gene sets
# Interesting/useful results: C2/CP:WIKIPATHWAYS,CP:KEGG,CP:PID, H
msigdbH <- msigdbr(species = "human", category = "C2", subcategory = "CP:PID")
# fixing format to work with fgsea
pathwaysH = split(x = msigdbH$gene_symbol, f = msigdbH$gs_name)

#DE_up <- subset(markers, avg_log2FC > 0)
#ranks <- DE_up$avg_log2FC
ranks <- nonT$avg_log2FC
names(ranks) <- rownames(nonT)
ranks <- sort(ranks, decreasing = T)
outgsea <- data.frame()
for (i in 1:length(pathwaysH)){
  outgsea <- rbind(outgsea,
                   fgsea(pathways = pathwaysH[i],
                         stats = ranks,
                         minSize=10,
                         maxSize=500,
                         nproc=1, # for some reason get an error if allowed to use more than one processor
                         nPermSimple = 5000))
}
head(outgsea[order(outgsea$padj),])



FeaturePlot(mSeurat, c("CD3D", "PITPNC1", "CAMK1D", "ATXN1"))
VlnPlot(mSeurat, "PITPNC1", split.by = "Tissue", split.plot = T)
VlnPlot(mSeurat, "ZSWIM6", split.by = "Group", split.plot = T)
# ZSWIM is in glial and neuronal cells -_- !!!
# PTPRC/CD45 is mature immune cells
FeaturePlot(mSeurat, c("CD74", "ZSWIM6"))
FeatureScatter(subset(mSeurat, seurat_clusters == 2), "ZSWIM6", "CD3D")



# exploratory plots {.tabset}

## plots

mSeurat$tissuegroup <- paste(mSeurat$Tissue, mSeurat$Group)
DimPlot(mSeurat, group.by = "tissuegroup", shuffle = T, cols = c("darkred", "orange", "steelblue1","purple"))


# View(as.matrix(mSeurat@assays$Protein@data))
# FeaturePlot(mSeurat, rownames(mSeurat@assays$Protein), slot = "data")
FeaturePlot(mSeurat, "PD1.P", slot = "data")
RidgePlot(mSeurat, "CD122.P")


FeaturePlot(mSeurat, "CXCR5")
FeaturePlot(mSeurat, "FCGR3A")
DimPlot(mSeurat, group.by = "n_chains")

RidgePlot(mSeurat, "Hashtag1")
RidgePlot(mSeurat, "Hashtag2")
RidgePlot(mSeurat, "Hashtag3")
RidgePlot(mSeurat, "Hashtag4")


## every protein

for (x in rownames(mSeurat@assays$Protein)){
  print(RidgePlot(mSeurat, x, log = T))
}


# barplots of composition {.tabset}

## clusters / donor


tissue_group <- unique(paste(mSeurat$Tissue, mSeurat$Group))


# do first for all cells (mSeurat)
for (y in tissue_group){
  df <- subset(mSeurat, 
               Tissue == strsplit(y, " ")[[1]][1] &
                 Group == strsplit(y, " ")[[1]][2])@meta.data
  output <- data.frame()
  
  # do the calculations for how many cells from each cluster in each donor
  for (x in unique(df$Patient)){
    input <- subset(df, Patient == x)
    input <- rename(count(input, clusters_named), Freq = n)
    input$Donor <- x
    output <- rbind(output, input)
  }
  
  # make the plot
  print(ggplot(output, aes(x = Donor, y = Freq, fill = clusters_named)) +
          geom_bar(position="fill", stat="identity", width = 0.5) +
          scale_fill_brewer(palette = "Set1") +
          theme_classic() +
          labs(subtitle = y) +
          theme(aspect.ratio = 1))
}



# now repeat for only cells with TCR (hasTCR)
for (y in tissue_group){
  df <- subset(hasTCR, 
               Tissue == strsplit(y, " ")[[1]][1] &
                 Group == strsplit(y, " ")[[1]][2])@meta.data
  output <- data.frame()
  
  # do the calculations for how many cells from each cluster in each donor
  for (x in unique(df$Patient)){
    input <- subset(df, Patient == x)
    input <- rename(count(input, clusters_named), Freq = n)
    input$Donor <- x
    output <- rbind(output, input)
  }
  
  # make the plot
  print(ggplot(output, aes(x = Donor, y = Freq, fill = clusters_named)) +
          geom_bar(position="fill", stat="identity", width = 0.5) +
          scale_fill_brewer(palette = "Set1") +
          theme_classic() +
          labs(subtitle = paste(y, "excluding cells with no TCR")) +
          theme(aspect.ratio = 1))
}

remove(x, y, tissue_group, df, output, input)


## clusters / pooled donors



tissue_group <- unique(paste(mSeurat$Tissue, mSeurat$Group))

output <- data.frame()
# do first for all cells (mSeurat)
for (y in tissue_group){
  df <- subset(mSeurat, 
               Tissue == strsplit(y, " ")[[1]][1] &
                 Group == strsplit(y, " ")[[1]][2])@meta.data
  
  df <- rename(count(df, seurat_clusters), Freq = n)
  df$tissue_group <- y
  output <- rbind(output, df)
}

# make the plot
print(ggplot(output, aes(x = tissue_group, y = Freq, fill = seurat_clusters)) +
        geom_bar(position="fill", stat="identity", width = 0.5) +
        scale_fill_brewer(palette = "Set1") +
        scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        labs(subtitle = "all cells"))



output <- data.frame()
# now repeat for only cells with TCR (hasTCR)
for (y in tissue_group){
  df <- subset(hasTCR, 
               Tissue == strsplit(y, " ")[[1]][1] &
                 Group == strsplit(y, " ")[[1]][2])@meta.data
  
  df <- rename(count(df, seurat_clusters), Freq = n)
  df$tissue_group <- y
  output <- rbind(output, df)
}

# make the plot
print(ggplot(output, aes(x = tissue_group, y = Freq, fill = seurat_clusters)) +
        geom_bar(position="fill", stat="identity", width = 0.5) +
        scale_fill_brewer(palette = "Set1") +
        scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        labs(subtitle = "only cells with TCR"))


remove(y, tissue_group, df, output)

## tissuegroup / cluster


output <- data.frame()
# do first for all cells (mSeurat)
for (y in unique(mSeurat$tissuegroup)){
  df <- subset(mSeurat, 
               Tissue == strsplit(y, " ")[[1]][1] &
                 Group == strsplit(y, " ")[[1]][2])@meta.data
  
  df <- rename(count(df, clusters_named), Freq = n)
  df$tissue_group <- y
  output <- rbind(output, df)
}

# make the plot
print(ggplot(output, aes(fill = tissue_group, y = Freq, x = clusters_named)) +
        geom_bar(position="fill", stat="identity", width = 0.5) +
        scale_fill_manual(values = c("darkred", "red", "steelblue1","lightblue")) +
        scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        RotatedAxis())

remove(output, y, df)




# diff expression between groups {.tabset}

## compare total CXCR5+/- in blood

# test for all cells
Idents(mSeurat) <- "Group"
diffgenes <- FindMarkers(subset(mSeurat, Tissue == "Blood"), 
                         ident.1 = "CXCR5positive", 
                         ident.2 = "CXCR5negative")
EnhancedVolcano(diffgenes, lab = paste0("italic('",rownames(diffgenes), "')"),
                parseLabels = T, 
                xlab = expression('Log'[2]*'(fold change)'),
                ylab = expression('-Log'[10]*'(adj. p-value)'),
                x = 'avg_log2FC', y = 'p_val_adj',
                pCutoff = 5e-2,
                FCcutoff = 0.25,
                drawConnectors = T,
                lengthConnectors = unit(0.01, "npc"),
                widthConnectors = unit(0.5, "npc"),
                colAlpha = 1,
                legendPosition = "none", caption = "all cells")

# test for without TCRless cells
Idents(hasTCR) <- "Group"
diffgenes <- FindMarkers(subset(hasTCR, Tissue == "Blood"), 
                         ident.1 = "CXCR5positive", 
                         ident.2 = "CXCR5negative")
EnhancedVolcano(diffgenes, lab = paste0("italic('",rownames(diffgenes), "')"),
                parseLabels = T, 
                xlab = expression('Log'[2]*'(fold change)'),
                ylab = expression('-Log'[10]*'(adj. p-value)'),
                x = 'avg_log2FC', y = 'p_val_adj',
                pCutoff = 5e-2,
                FCcutoff = 0.25,
                drawConnectors = T,
                lengthConnectors = unit(0.01, "npc"),
                widthConnectors = unit(0.5, "npc"),
                colAlpha = 1,
                legendPosition = "none", caption = "only cells with TCR")



## compare total CXCR5+/- in tonsil

# test for all cells
Idents(mSeurat) <- "Group"
diffgenes <- FindMarkers(subset(mSeurat, Tissue == "Tonsil"), 
                         ident.1 = "CXCR5positive", 
                         ident.2 = "CXCR5negative")
EnhancedVolcano(diffgenes, lab = paste0("italic('",rownames(diffgenes), "')"),
                parseLabels = T, 
                xlab = expression('Log'[2]*'(fold change)'),
                ylab = expression('-Log'[10]*'(adj. p-value)'),
                x = 'avg_log2FC', y = 'p_val_adj',
                pCutoff = 5e-2,
                FCcutoff = 0.25,
                drawConnectors = T,
                lengthConnectors = unit(0.01, "npc"),
                widthConnectors = unit(0.5, "npc"),
                colAlpha = 1,
                legendPosition = "none", caption = "all cells")

# test for without TCRless cells
Idents(hasTCR) <- "Group"
diffgenes <- FindMarkers(subset(hasTCR, Tissue == "Tonsil"), 
                         ident.1 = "CXCR5positive", 
                         ident.2 = "CXCR5negative")
EnhancedVolcano(diffgenes, lab = paste0("italic('",rownames(diffgenes), "')"),
                parseLabels = T, 
                xlab = expression('Log'[2]*'(fold change)'),
                ylab = expression('-Log'[10]*'(adj. p-value)'),
                x = 'avg_log2FC', y = 'p_val_adj',
                pCutoff = 5e-2,
                FCcutoff = 0.25,
                drawConnectors = T,
                lengthConnectors = unit(0.01, "npc"),
                widthConnectors = unit(0.5, "npc"),
                colAlpha = 1,
                legendPosition = "none", caption = "only cells with TCR")
```


## compare total CXCR5+/- protein in blood

# test for all cells
Idents(mSeurat) <- "Group"
diffgenes <- FindMarkers(subset(mSeurat, Tissue == "Blood"), 
                         ident.1 = "CXCR5positive", 
                         ident.2 = "CXCR5negative", assay = "Protein")
RidgePlot(mSeurat, rownames(diffgenes), log = T) &
  theme(text = element_text(size=6))


# test for cells having TCR
Idents(hasTCR) <- "Group"
diffgenes <- FindMarkers(subset(hasTCR, Tissue == "Blood"), 
                         ident.1 = "CXCR5positive", 
                         ident.2 = "CXCR5negative", assay = "Protein")
RidgePlot(hasTCR, rownames(diffgenes), log = T) &
  theme(text = element_text(size=6))



## compare total CXCR5+/- protein in tonsil

# test for all cells
Idents(mSeurat) <- "Group"
diffgenes <- FindMarkers(subset(mSeurat, Tissue == "Tonsil"), 
                         ident.1 = "CXCR5positive", 
                         ident.2 = "CXCR5negative", assay = "Protein")
RidgePlot(mSeurat, rownames(diffgenes), log = T, assay = "Protein") &
  theme(text = element_text(size=6))


# test for cells having TCR
Idents(hasTCR) <- "Group"
diffgenes <- FindMarkers(subset(hasTCR, Tissue == "Tonsil"), 
                         ident.1 = "CXCR5positive", 
                         ident.2 = "CXCR5negative", assay = "Protein")
RidgePlot(mSeurat, rownames(diffgenes), log = T, assay = "Protein") &
  theme(text = element_text(size=6))


# PD-1 expression {.tabset}

## in general


# comparison of PDCD1 transcript expression
VlnPlot(mSeurat, "PDCD1")

# show PD-1 expression via protein
RidgePlot(mSeurat, "PD1.P", log = T, group.by = "clusters_named") + 
  geom_vline(xintercept = 1.6) +
  NoLegend()

# setup variables
marker <- c("PDCD1", "PD1.P")
threshold <- c(0, 1.6)

# pie charts of proportion of cluster with PDCD1 expression
for (i in 1:length(marker)){
  is_positive <- as.vector(
    ifelse(FetchData(mSeurat, vars = marker[i]) > threshold[i], "Positive", "Negative"))
  df <- data.frame("marker" = is_positive, "cluster" = mSeurat$clusters_named)
  df <- rename(count(df, cluster, marker), Freq = n)
  plots <- list()
  for (x in unique(df$cluster)){
    input <- subset(df, cluster == x)
    plots[[x]] <- ggplot(input, aes(x="", y=Freq, fill=marker)) +
      geom_bar(stat="identity", width=1) +
      coord_polar("y", start=0) +
      theme_void() +
      labs(subtitle = paste(x, ":", marker[i])) +
      theme(text = element_text(size=6))
  }
  do.call("grid.arrange", c(plots, ncol=3))
}

remove(marker, threshold, i, df, x, input, plots, is_positive)
```


## within CXCR5+/-
```{r PD1 and CXCR5, warning=FALSE, message=FALSE}

VlnPlot(mSeurat, "PDCD1", group.by = "Tissue", split.by = "Group", split.plot = T)
VlnPlot(mSeurat, "CXCR5", group.by = "Tissue", split.by = "Group", split.plot = T)

group <- c("CXCR5positive", "CXCR5negative")
tissue <- c("Blood", "Tonsil")
marker <- c("PDCD1", "CXCR5", "PD1.P")
threshold <- c(0,0, 1.6)


for (i in 1:length(marker)){
  plots <- list()
  for (x in group){
    for (y in tissue){
      data <- subset(mSeurat, Group == x & Tissue == y)
      is_positive <- as.vector(
        ifelse(FetchData(data, vars = marker[i]) > threshold[i], "Positive", "Negative"))
      df <- data.frame("marker" = is_positive, sample = paste(y, x))
      df <- rename(count(df, marker, sample), Freq = n)
      plots[[paste(y,x)]] <- ggplot(df, aes(x="", y=Freq, fill=marker)) +
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) +
        theme_void() +
        labs(subtitle = paste(y, x, ":", marker[i]))
    }
  }
  do.call("grid.arrange", c(plots, ncol=2))
}

remove(group, tissue, marker, threshold, i, x, y, data, is_positive, df, plots)

