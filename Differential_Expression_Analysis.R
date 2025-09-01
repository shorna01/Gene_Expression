# Step 2: Differential Expression Analysis


if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
library(DESeq2)


# Re-create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = round(expr_filtered), # must be integers
                              colData = metadata,
                              design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)

library(ashr)

# Extract results (Long_COVID vs Control)
res <- results(dds, contrast = c("condition", "Long_COVID", "Control"))
res <- lfcShrink(dds, contrast = c("condition", "Long_COVID", "Control"),
                 res = res, type = "ashr")   # shrink fold-changes

# Order by adjusted p-value
res <- res[order(res$padj), ]
head(res)


#Installing EnhancedVolcano Plot

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)


# Volcano plot

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.001,
                FCcutoff = 2,
                title = "Volcano_Plot of DEG's")



# Heatmap (Top 20 DEGs)

library(pheatmap)

topgenes <- rownames(res)[1:20]
mat <- assay(vst(dds))[topgenes, ]   # variance stabilized counts
pheatmap(mat, annotation_col = metadata, show_rownames = TRUE,
         main = "Top 20 Different Expressed Genes Heatmap")


#  Results

write.csv(as.data.frame(res), file = "DEGs_full.csv")

# Convert DESeqResults to data frame and remove NAs
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj), ]

# Top 15 up-regulated
top_up <- head(res_df[order(res_df$log2FoldChange, decreasing = TRUE), ], 15)
write.csv(top_up, file = "Top15_Up regulated.csv", row.names = TRUE)

# Top 15 down-regulated
top_down <- head(res_df[order(res_df$log2FoldChange, decreasing = FALSE), ], 15)
write.csv(top_down, file = "Top15_Down regulated.csv", row.names = TRUE)


#Location of saved 15 files
getwd()
setwd("C:/Users/DELL/Desktop/Gene Expression Biomarkers/Project_New")



#Checking 15 files(Top & Bottom)

file.exists("Top15_Down regulated.csv")
file.exists("Top15_Up regulated.csv")







######## Protein-Protein_Interaction_Network_Analysis ######


# Significant DEGs (padj < 0.05 and |log2FC| > 1)
sig_DEGs <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sig_DEGs), "Significant_DEGs.csv")


# Read hub genes exported from Cytoscape (if saved as CSV)
hub_genes <- c("LHX5","IGF2","SLC22A18AS","C16orf71",
               "PBX4","PKNOX2","FEZF2","RELN","MEIS2","PBX2")

# Subset your DESeq2 results for these hub genes
hub_res <- res_df[rownames(res_df) %in% hub_genes, ]

# Check the table
hub_res

# Plot log2 fold change
library(ggplot2)
ggplot(hub_res, aes(x=reorder(rownames(hub_res), log2FoldChange),
                    y=log2FoldChange)) +
  geom_bar(stat="identity", fill="green") +
  coord_flip() +
  xlab("Hub Genes") +
  ylab("log2 Fold Change") +
  ggtitle("DEGs Hub Genes log2FC")












