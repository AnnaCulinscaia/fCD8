
# Libraries ---------------------------------------------------------------

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
                x = 'avg_log2FC', y = 'p_val_adj',
                pCutoff = 5e-2,
                FCcutoff = 0.25,
                drawConnectors = T,
                lengthConnectors = unit(0.01, "npc"),
                widthConnectors = unit(0.5, "npc"),
                colAlpha = 1,
                legendPosition = "none", caption = "all cells")
ggsave(file = "plotplease.pdf", plot = plot1, device = NULL, width = 20, height = 20, units = "cm")
