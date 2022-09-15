
# Libraries ---------------------------------------------------------------
library(DESeq2)
library(EnhancedVolcano)
#•	Do a differential expression between all CXCR5+ and CXCR5- in blood
Idents(mSeurat) <- "Group"
diffgenes <- FindMarkers(subset(mSeurat, Tissue == "Blood"), 
                         ident.1 = "CXCR5positive", 
                         ident.2 = "CXCR5negative")

plot1<- EnhancedVolcano(diffgenes, lab = paste0("italic('",rownames(diffgenes), "')"),
                parseLabels = T, 
                xlab = expression('Log'[2]*'(fold change)'),
                ylab = expression('-Log'[10]*'(adj. p-value)'),
                title = 'CXCR5+ vs CXCR5- in blood',
                x = 'avg_log2FC', y = 'p_val_adj',
                pCutoff = 5e-2,
                FCcutoff = 0.25,
                drawConnectors = T,
                lengthConnectors = unit(0.01, "npc"),
                widthConnectors = unit(0.5, "npc"),
                colAlpha = 1,
                legendPosition = "none", caption = "all cells")
ggsave(file = "Volcano_CXCR5+_vs_CXCR5_blood.png", plot = plot1, device = NULL, width = 25, height = 13, units = "cm")

#differential expression between all CXCR5+ and CXCR5- in tonsil
Idents(mSeurat) <- "Group"
diffgenes2 <- FindMarkers(subset(mSeurat, Tissue == "Tonsil"), 
                         ident.1 = "CXCR5positive", 
                         ident.2 = "CXCR5negative")

plot2<- EnhancedVolcano(diffgenes2, lab = paste0("italic('",rownames(diffgenes2), "')"),
                        parseLabels = T, 
                        title = 'CXCR5+ vs CXCR5- in tonsil',
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
ggsave(file = "Volcano_CXCR5+_vs_CXCR5_tonsil.png", plot = plot2, device = NULL, width = 25, height = 13, units = "cm")

#differential expression between CXCR5+ and CXCR5- in blood per donor
make_volcano <- function(dnr, tis){
  sSeurat <- subset(mSeurat, Tissue == tis & Patient == dnr)
  diffgenes<- FindMarkers(sSeurat, ident.1 = "CXCR5positive", ident.2 = "CXCR5negative")
  plot<- EnhancedVolcano(diffgenes, lab = paste0("italic('",rownames(diffgenes), "')"),
                          parseLabels = T, 
                          title = paste('CXCR5+ vs CXCR5- for',dnr,'`s',tis),
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
  ggsave(file = paste0("Volcano_CXCR5+_vs_CXCR5_",dnr,"_", tis, ".png", collapse =NULL), plot = plot, device = NULL, width = 25, height = 13, units = "cm")
  
}

donor_list <- c("Dnr58","Dnr69","Dnr72","Dnr83" )
tissue_list <- c("Blood", "Tonsil")

for (tis_type in tissue_list){
  for (dnr_name in donor_list){
    make_volcano(dnr_name, tis_type)
  }
}
#•	Do a differential expression between CXCR5+ in blood and tonsil
Idents(mSeurat) <- "Tissue"

diffgenes3 <- FindMarkers(subset(mSeurat, Group == "CXCR5positive"), 
                          ident.1 = "Blood", 
                          ident.2 = "Tonsil")

plot3<- EnhancedVolcano(diffgenes3, lab = paste0("italic('",rownames(diffgenes3), "')"),
                        parseLabels = T, 
                        title = 'Blood vs Tonsil in CXCR5+',
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
ggsave(file = "Volcano_Blood_vs_Tonsil_CXCR5+.png", plot = plot3, device = NULL, width = 25, height = 13, units = "cm")
#•	Do a differential expression between CXCR5- in blood and tonsil
Idents(mSeurat) <- "Tissue"

diffgenes4 <- FindMarkers(subset(mSeurat, Group == "CXCR5negative"), 
                          ident.1 = "Blood", 
                          ident.2 = "Tonsil")

plot4<- EnhancedVolcano(diffgenes4, lab = paste0("italic('",rownames(diffgenes4), "')"),
                        parseLabels = T, 
                        title = 'Blood vs Tonsil in CXCR5-',
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
ggsave(file = "Volcano_Blood_vs_Tonsil_CXCR5-.png", plot = plot4, device = NULL, width = 25, height = 13, units = "cm")
