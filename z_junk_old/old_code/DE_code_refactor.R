# ====================================================================
# Initialisation & Libraries
# ====================================================================

library(DESeq2)
library(ggplot2)
library(WGCNA)
library(pheatmap)
library(vegan) 
library(pairwiseAdonis)
library(glue)
library(dplyr)
library(purrr)
library(readr)
library(stringr) 
library(tibble) 
library(apeglm)
library(gtools)
library(clusterProfiler)
library(org.Mm.eg.db)  
library(enrichplot)

# ====================================================================
# Phase 1: Preparing Your Data & DESeq2 Setup
# ====================================================================

# 1. Load counts using the Jupyter notebook file structure
count_files <- gtools::mixedsort(Sys.glob("o_outputs/sample_*/gene_counts.txt"))

# Extract the Geneid column from the first file to initialise the master matrix
master_counts <- read.table(count_files[1], header = TRUE, skip = 1, stringsAsFactors = FALSE)[, 1, drop = FALSE]
colnames(master_counts) <- "Geneid"

# Loop through and append count columns, dynamically naming them by folder
for (file in count_files) {
    sample_id <- basename(dirname(file))
    master_counts[[sample_id]] <- read.table(file, header = TRUE, skip = 1)[, 7]
}

# 2. Clean up the featureCounts meta-rows at the bottom (starting with double underscores)
master_counts <- master_counts %>% filter(!str_detect(Geneid, "^__"))

# 3. Set Geneid as the row names (required by DESeq2)
counts_final <- master_counts
rownames(counts_final) <- counts_final$Geneid
counts_final$Geneid <- NULL

# 4. Check matrix dimensions and save
# print(dim(counts_final)) 
write.csv(counts_final, file = "master_count_matrix.csv", row.names = TRUE)

# ====================================================================
# Phase 2: Variance Stabilizing Transformation (VST) and PCA plot
# ====================================================================

# 1. Dynamically extract the sample number and assign the correct group
sample_numbers <- as.numeric(str_extract(colnames(counts_final), "\\d+"))

sample_info <- data.frame(
  condition = case_when(
    sample_numbers %in% 1:3   ~ "pulldown_unlab",
    sample_numbers %in% 4:6   ~ "pulldown_lab",
    sample_numbers %in% 7:9   ~ "input_unlab",
    sample_numbers %in% 10:12 ~ "input_lab",
    TRUE                      ~ "unknown_sample"
  ),
  row.names = colnames(counts_final)
)

# 2. Convert to factor and explicitly set the baseline reference level
sample_info$condition <- factor(
  sample_info$condition, 
  levels = c("input_unlab", "input_lab", "pulldown_unlab", "pulldown_lab")
)

# 3. Build the DESeq2 Dataset Object
dds <- DESeqDataSetFromMatrix(
  countData = counts_final,
  colData = sample_info,
  design = ~ condition
)

# 4. Broad Filter: Keep genes with >= 10 reads in at least 3 of the 12 samples
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

# 5. Apply Variance Stabilizing Transformation (VST)
vsd <- vst(dds, blind = TRUE)

# 6. Generate the PCA Data Structure
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# # 7. Plot using ggplot2
# ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
#   geom_point(size = 5, alpha = 0.8) +  
#   xlab(paste0("PC1: ", percent_var[1], "% variance")) +
#   ylab(paste0("PC2: ", percent_var[2], "% variance")) +
#   coord_fixed() +
#   theme_minimal() +
#   scale_color_brewer(palette = "Set1") + 
#   labs(
#     title = "RNA-Seq Raw Data PCA Plot (VST Normalised)",
#     subtitle = "Comparing labelled vs. unlabelled samples across input and pulldown fractions",
#     color = "Biological Group"
#   ) +
#   theme(
#     plot.title = element_text(face = "bold", size = 14),
#     axis.title = element_text(size = 12),
#     legend.text = element_text(size = 10)
#   ) 

# ggsave(
#   filename = "PCA_plot_12samples.png", 
#   plot = last_plot(),                  
#   device = "png", 
#   width = 8, 
#   height = 6, 
#   dpi = 300                            
# )

# cat("PCA plot successfully saved to your directory!\n") 

# ====================================================================
# Phase 3: Differential Expression Analysis (DESeq2)
# ====================================================================

# 1. Explicitly set 'pulldown_unlab' as the reference level
dds$condition <- relevel(dds$condition, ref = "pulldown_unlab")

# 2. Run the core pipeline 
dds <- DESeq(dds)
print(resultsNames(dds))

# 3. Extract the specific comparison
res_enrichment <- results(
  dds, 
  contrast = c("condition", "pulldown_lab", "pulldown_unlab"),
  alpha = 0.05
)

# 4. Apply apeglm log2 fold change shrinkage 
res_shrink <- lfcShrink(
  dds, 
  coef = "condition_pulldown_lab_vs_pulldown_unlab", 
  type = "apeglm"
)

summary(res_enrichment)

# Extract results with a relaxed independent filtering threshold
res_relaxed <- results(dds, contrast = c("condition", "pulldown_lab", "pulldown_unlab"), alpha = 0.1)

# Apply shrinkage to this new result
res_relaxed_shrunk <- lfcShrink(dds, coef = "condition_pulldown_lab_vs_pulldown_unlab", type = "apeglm")

# Filter for targets that are enriched (LFC > 0) at padj < 0.1
expanded_hits <- as.data.frame(res_relaxed_shrunk) %>%
  rownames_to_column(var = "gene_id") %>%
  filter(padj < 0.1 & log2FoldChange > 0)

message("Number of hits at padj < 0.1: ", nrow(expanded_hits)) 

# Convert your shrunk results into a standard data frame
res_df <- as.data.frame(res_shrink) %>%
  rownames_to_column(var = "gene_id")

# ==============================================================================
# GENERATE TIER 1: Core High-Confidence Hits 
# ==============================================================================
# Criteria: Statistically significant after multi-test correction, and enriched
tier1_core_hits <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  arrange(padj) 

write.csv(tier1_core_hits, file = "pulldown_Tier1_core_hits_padj05.csv", row.names = FALSE)

# ==============================================================================
# GENERATE TIER 2: Broad Candidate Hits 
# ==============================================================================
# Criteria: Significant by raw p-value, enriched, and excludes Tier 1 to avoid duplication
tier2_candidate_hits <- res_df %>%
  filter(pvalue < 0.05 & log2FoldChange > 0) %>%
  filter(!gene_id %in% tier1_core_hits$gene_id) %>% 
  arrange(pvalue) 

write.csv(tier2_candidate_hits, file = "pulldown_Tier2_candidate_hits_pvalue05.csv", row.names = FALSE)

message("Files successfully written!")
message("Tier 1 (Strict padj < 0.05): ", nrow(tier1_core_hits), " genes.")
message("Tier 2 (Raw pvalue < 0.05, excluding Tier 1): ", nrow(tier2_candidate_hits), " genes.") 

# ==============================================================================
# GENE ONTOLOGY (GO) ENRICHMENT
# ==============================================================================

# 1. Prepare your combined gene list
target_genes <- c(tier1_core_hits$gene_id, tier2_candidate_hits$gene_id)

# 2. Run the Gene Ontology Enrichment for Mouse
ego_mouse <- enrichGO(gene          = target_genes,
                      OrgDb         = org.Mm.eg.db,  
                      keyType       = "SYMBOL",      
                      ont           = "BP",          
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)

# 3. Save the mouse results table to a CSV file
go_mouse_results <- as.data.frame(ego_mouse)
write.csv(go_mouse_results, file = "mouse_pulldown_GO_results.csv", row.names = FALSE)

# ==============================================================================
# VISUALISE THE MOUSE HITS
# ==============================================================================

# 1. Dotplot (Shows ratio of genes vs statistical significance)
print(dotplot(ego_mouse, showCategory = 15) + ggtitle("Top 15 Mouse Biological Processes"))

# 2. Barplot (Shows absolute gene counts)
print(barplot(ego_mouse, showCategory = 15) + ggtitle("Top 15 Mouse Biological Processes"))

# ==============================================================================
# CILIA SPECIFIC FILTERING
# ==============================================================================

go_results <- as.data.frame(ego_mouse)

# Filter for the specific ciliary terms
cilia_terms <- go_results %>%
  filter(grepl("cilium assembly|cilium organization", Description, ignore.case = TRUE))

print(cilia_terms[, c("ID", "Description", "Count")])

# Extract and clean the gene list
cilia_genes_list <- cilia_terms$geneID %>%
  paste(collapse = "/") %>%       
  strsplit(split = "/") %>%       
  unlist() %>%                    
  unique() %>%                    
  sort()                          

message("\nFound ", length(cilia_genes_list), " unique genes matching cilia terms:")
print(cilia_genes_list)