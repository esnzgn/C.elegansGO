# Load necessary libraries
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(dplyr)
library(ReactomePA)
library(org.Ce.eg.db)


# worm list enrichment
# reading
file_path <- "./06-2024//List of verified genes.xlsx"
gene_df <- read_excel(file_path)
head(gene_df)  # Inspect the structure
# bp enrichment
symbol_col <- grep("WormBase Gene", colnames(gene_df), ignore.case = TRUE, value = TRUE)

gene_df <- gene_df %>%
  rename(GeneID = all_of(symbol_col))

# Step 2: Convert WormBase IDs to ENTREZID
gene_map <- bitr(gene_df$GeneID,
                 fromType = "WORMBASE",    # WormBase Gene ID
                 toType = "ENTREZID",
                 OrgDb = org.Ce.eg.db)

# entrez_ids <- unique(gene_map$ENTREZID)
#
# # Step 3: GO enrichment (M F)
# ego_ce_mf <- enrichGO(gene          = entrez_ids,
#                       OrgDb         = org.Ce.eg.db,
#                       keyType       = "ENTREZID",
#                       ont           = "MF",
#                       pAdjustMethod = "BH",
#                       pvalueCutoff  = 0.05,
#                       qvalueCutoff  = 0.2,
#                       readable      = TRUE)

# head(ego_ce_mf)
# # Step 4: Plot results
# dotplot(ego_ce_mf, showCategory = 20) + ggtitle("GO Enrichment (C. elegans - MF)")

# Step 3: GO enrichment (Biological Process)
ego_ce_bp <- enrichGO(gene          = entrez_ids,
                      OrgDb         = org.Ce.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE)

# Step 4: Plot results
dotplot(ego_ce_bp, showCategory = 20) + ggtitle("GO Enrichment (C. elegans - BP)")

# KEGG uses Entrez IDs, so you may need to convert:
ekegg_ce <- enrichKEGG(gene         = entrez_ids,
                       organism     = 'cel',
                       pvalueCutoff = 0.05)

# dotplot(ekegg_ce, showCategory = 20) + ggtitle("KEGG Enrichment (C. elegans)")

# write.csv(as.data.frame(ego_ce_bp), "Celegans_GO_BP.csv", row.names = FALSE)
# write.csv(as.data.frame(ekegg_ce), "Celegans_KEGG.csv", row.names = FALSE)

# Step 1: Load the gene list
file_path <- "./30-05-25/List of human orthologs_RNAi screen.xlsx"
gene_df <- read_excel(file_path)
head(gene_df)  # Inspect the structure

# Step 2: Extract and convert gene symbols to Entrez IDs
# Replace 'GeneSymbolColumn' with the actual column name from the file
symbol_col <- grep("symbol", names(gene_df), ignore.case = T, value = T)

gene_df <- gene_df %>%
  rename(GeneSymbolColumn = all_of(symbol_col))

gene_symbols <- gene_df$GeneSymbolColumn  # Update with the correct column name
entrez_ids <- bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez_ids <- unique(entrez_ids$ENTREZID)

# Step 3: Perform GO enrichment (Biological Process)
ego_bp <- enrichGO(gene         = entrez_ids,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "ENTREZID",
                   ont          = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable     = TRUE)

ego_mf <- enrichGO(gene         = entrez_ids,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "ENTREZID",
                   ont          = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable     = TRUE)

# Step 4: Perform KEGG pathway enrichment
ekegg <- enrichKEGG(gene         = entrez_ids,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)

# Step 5: Visualization
head(ego_bp)
dotplot(ego_bp, showCategory = 20) + ggtitle("GO Enrichment (Biological Process)")
head(ego_mf)
dotplot(ego_mf, showCategory = 20) + ggtitle("GO Enrichment (Moleculal Function)")
# head(ekegg)
# dotplot(ekegg, showCategory = 20) + ggtitle("KEGG Pathway Enrichment")

# reactome instead of kegg
ereactome <- enrichPathway(gene = entrez_ids,
                           organism = "human",
                           pvalueCutoff = 0.05,
                           readable = TRUE)
head(ereactome)
dotplot(ereactome, showCategory = 20) + ggtitle("Reactome Pathway Enrichment")


# Optional: Save results
# write.csv(as.data.frame(ego_bp), "GO_BP_Enrichment.csv", row.names = FALSE)
# write.csv(as.data.frame(ekegg), "KEGG_Enrichment.csv", row.names = FALSE)
