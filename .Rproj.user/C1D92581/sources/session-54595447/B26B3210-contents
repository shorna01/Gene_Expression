########Step-4#########

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")



# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)



#(A) Pathway Enrichment for All Significant DEGs#

# Convert gene symbols into Entrez IDs

sig_gene_symbols <- rownames(sig_DEGs)
sig_gene_entrez <- bitr(sig_gene_symbols,
                        fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)

#  Gene Ontology Enrichment
ego_all <- enrichGO(gene          = sig_gene_entrez$ENTREZID,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "ALL",   # BP, MF, CC, ALL
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)

dotplot(ego_all, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") +
  ggtitle(" Gene Ontology Enrichment (All DEGs)")

# Save results
write.csv(as.data.frame(ego_all), "GO_All_DEGs.csv")


# EGG Enrichment
ekegg_all <- enrichKEGG(
  gene          = sig_gene_entrez$ENTREZID,
  organism      = "hsa",
  pvalueCutoff  = 0.1,   # relaxed cutoff
  qvalueCutoff  = 0.2
)

if (nrow(as.data.frame(ekegg_all)) > 0) {
  barplot(ekegg_all, showCategory=20, title="KEGG Pathways (All DEGs)")
  dotplot(ekegg_all, showCategory=20, title="KEGG Pathways (All DEGs)")
  write.csv(as.data.frame(head(ekegg_all, 10)),
            "KEGG_Top10.csv", row.names = FALSE)
} else {
  message("⚠ No KEGG enrichment found.")
}


# eactome Enrichment
ereact_all <- enrichPathway(
  gene          = sig_gene_entrez$ENTREZID,
  organism      = "human",
  pvalueCutoff  = 0.1,
  readable      = TRUE
)

if (nrow(as.data.frame(ereact_all)) > 0) {
  dotplot(ereact_all, showCategory=20, title="Reactome Pathways (All DEGs)")
  barplot(ereact_all, showCategory=20, title="Reactome Pathways (All DEGs)")
  write.csv(as.data.frame(head(ereact_all, 10)),
            "Reactome_Top10.csv", row.names = FALSE)
} else {
  message("⚠ No Reactome enrichment found.")
}



####(B) Pathway Enrichment for Hub Genes Only

# Convert hub genes → Entrez IDs
hub_gene_entrez <- bitr(hub_genes,
                        fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)

# Gene Ontology Enrichment
ego_hub <- enrichGO(gene          = hub_gene_entrez$ENTREZID,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)

dotplot(ego_hub, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") +
  ggtitle("Gene Ontology Enrichment (Hub Genes)")

# Save results
write.csv(as.data.frame(ego_hub), "Gene_Ontology _HubGenes.csv")

# EGG Enrichment
ekegg_hub <- enrichKEGG(gene          = hub_gene_entrez$ENTREZID,
                        organism      = 'hsa',
                        pvalueCutoff  = 0.05)
barplot(ekegg_hub, showCategory=15, title="KEGG Pathways (Hub Genes)")

# Save results
write.csv(as.data.frame(ekegg_hub), "KEGG_HubGenes.csv")

# ⚡ Reactome Enrichment
ereact_hub <- enrichPathway(gene         = hub_gene_entrez$ENTREZID,
                            organism     = "human",
                            pvalueCutoff = 0.05,
                            readable     = TRUE)
dotplot(ereact_hub, showCategory=15, title="Reactome Pathways (Hub Genes)")

# Save results
write.csv(as.data.frame(ereact_hub), "Reactome_HubGenes.csv")