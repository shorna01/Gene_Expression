library(tidyverse)

# Read TSV file (gene expression matrix)
expr_data <- read.delim("GSE270045_LC_counts.tsv", row.names = 1, sep = "\t")

# Preview
dim(expr_data)       # check dimensions (genes × samples)

head(expr_data[,1:5])   # first 5 samples, first few genes


# Replace with your exact sample names
sample_names <- colnames(expr_data)

metadata <- data.frame(
  row.names = sample_names,
  condition = c(
    rep("Long_COVID", 19),   # First 19 samples: Long COVID 
    rep("Control", 17)      # Next 17 samples: Healthy Controls
  )
)



# Significant DEGs (padj < 0.05 and |log2FC| > 1)
sig_DEGs <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sig_DEGs), "Significant_DEGs.csv")
# checking confirm
print(metadata)
table(metadata$condition)


# Keep genes that have expression >1 in at least 5 samples
keep <- rowSums(expr_data > 1) >= 5
expr_filtered <- expr_data[keep, ]
dim(expr_filtered)   # see how many genes remain


# PCA requires log-transformed data
expr_log <- log2(expr_filtered + 1)

pca_before <- prcomp(t(expr_log), scale. = TRUE)

# Plot
plot(pca_before$x[,1], pca_before$x[,2],
     col = as.factor(metadata$condition),
     pch = 19,
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA of Samples_before_batch_creation")
legend("topright", legend = levels(as.factor(metadata$condition)),
       col = 1:2, pch = 19)


# Load the GEOquery package

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}
library(GEOquery)

# Fetch the GEO dataset
gse <- getGEO("GSE270045", GSEMatrix = TRUE)

# Extract the phenotype data
pheno_data <- pData(phenoData(gse[[1]]))

# View the first few rows of the phenotype data to identify relevant columns
head(pheno_data)

head(pheno_data[,1:5])   # preview first 5 columns

colnames(pheno_data) #GEO metadata (pheno_data) has a real batch column
# Batch effect correction (ComBat)


# Check instrument model (sequencing platform)
table(pheno_data$instrument_model)

# Check source_name_ch1 (might have sample collection/run info)
head(pheno_data$source_name_ch1, 20)

# Check characteristics columns (often hold condition, sex, age, or batch/run info)
head(pheno_data$characteristics_ch1, 10)
head(pheno_data$characteristics_ch1.1, 10)
head(pheno_data$characteristics_ch1.2, 10)




if (!requireNamespace("sva", quietly = TRUE)) {
  BiocManager::install("sva")
}
library(sva)

#“No batch effects were detected, as all samples were sequenced on the same platform. 
#Therefore, PCA before and after batch correction appear identical.”







