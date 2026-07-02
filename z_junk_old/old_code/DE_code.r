# # Install required packages if missing ----
#  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "WGCNA", "pheatmap", "apeglm", "vegan", "pairwiseAdonis"))

# Input libraries -- assumes they are installed, if not use install.packages("NAME") or BiocManager::install("NAME")
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
library(httpgd) 
library(tibble) 
library(apeglm)

# ====================================================================
# Phase 1: Preparing Your Data & DESeq2 Setup
# ====================================================================

# 1. Point to the directory containing your individual featureCounts files
# Change this path to your actual counts folder
count_dir <- "count_data"

# 2. List all the count files (assuming they end in .txt)
count_files <- list.files(path = count_dir, pattern = "\\.txt$", full.names = TRUE)

# 3. Read and merge all files using an outer join (full_join)
master_matrix <- count_files %>%
  map(function(file) {
    # Extract a clean sample name from the filename to use as the column header
    sample_name <- basename(file) %>% 
      str_remove("_gene_counts.txt") %>%  # Clean up the suffix
      str_remove(".txt")                  # Catch any other text extensions
    
    # Read the file. featureCounts output has a 1-line header we need to skip,
    # and the columns are: Geneid, Chr, Start, End, Strand, Length, and the Count column.
    read_tsv(file, comment = "#", show_col_types = FALSE) %>%
      dplyr::select(Geneid, Count = last_col()) %>% # Grab the Gene ID and the very last column (the counts)
      rename(!!sample_name := Count)        # Dynamically rename the count column to the sample name
  }) %>%
  reduce(full_join, by = "Geneid")

# 4. Replace any NA values (genes missing in some samples) with 0
master_matrix[is.na(master_matrix)] <- 0

# 5. Clean up the featureCounts meta-rows at the bottom (they all start with double underscores)
# This removes __no_feature, __ambiguous, __alignment_not_unique, etc.
master_matrix <- master_matrix %>%
  filter(!str_detect(Geneid, "^__"))

# 6. Optional: Set Geneid as the row names (required by DESeq2)
counts_final <- as.data.frame(master_matrix)
rownames(counts_final) <- counts_final$Geneid
counts_final$Geneid <- NULL

# 7. Check matrix dimensions and preview the first few rows
head(counts_final)
dim(counts_final) # This will show you the exact same number of rows for every sample now

# 8. Save it for downstream differential expression
write.csv(counts_final, file = file.path(count_dir, "master_count_matrix.csv"), row.names = TRUE)

# ====================================================================
# Phase 2: Variance Stabilizing Transformation (VST) and PCA plot
# ====================================================================

# 1. Dynamically extract the sample number and assign the correct group
# This works completely independent of how R shuffles the row order!
sample_numbers <- as.numeric(str_extract(colnames(counts_final), "\\d+"))

sample_info <- data.frame(row.names = colnames(counts_final)) %>%
  mutate(
    condition = case_when(
      sample_numbers %in% 1:3   ~ "pulldown_unlab",
      sample_numbers %in% 4:6   ~ "pulldown_lab",
      sample_numbers %in% 7:9   ~ "input_unlab",
      sample_numbers %in% 10:12 ~ "input_lab",
      TRUE                      ~ "unknown_sample"
    )
)
  row.names = colnames(counts_final) # Binds the group labels directly to your sample columns


# 2. Convert to factor and explicitly set the baseline reference level
# Changing levels controls how DESeq2 evaluates comparisons down the road
sample_info$condition <- factor(
  sample_info$condition, 
  levels = c("input_unlab", "input_lab", "pulldown_unlab", "pulldown_lab")
)

# Quick verification check to make sure your columns matched up correctly
print(sample_info)

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

# 7. Plot using ggplot2
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 5, alpha = 0.8) +  # Made dots slightly larger for 12 samples
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  coord_fixed() +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") + # Beautiful, distinct color palette for 4 groups
  labs(
    title = "RNA-Seq Raw Data PCA Plot (VST Normalized)",
    subtitle = "Comparing labeled vs. unlabeled samples across input and pulldown fractions",
    color = "Biological Group"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) 

  ggsave(
  filename = "PCA_plot_12samples.png", 
  plot = last_plot(),                  # Grabs the ggplot you just made
  device = "png", 
  width = 8, 
  height = 6, 
  dpi = 300                            # High resolution for presentations/manuscripts
)

cat("Plot successfully saved to your directory!\n") 

# ====================================================================
# Phase 3: Differential Expression Analysis (DESeq2)
# ====================================================================

# 1. Run the core DESeq differential expression pipeline
# 1. Explicitly set 'pulldown_unlab' as the reference level
dds$condition <- relevel(dds$condition, ref = "pulldown_unlab")

# 2. Run the core pipeline 
# Now, the default comparison built into the model is lab vs unlab
dds <- DESeq(dds)

# Verify the name—it should now print "condition_pulldown_lab_vs_pulldown_unlab"
print(resultsNames(dds))

# 3. Extract the specific comparison
res_enrichment <- results(
  dds, 
  contrast = c("condition", "pulldown_lab", "pulldown_unlab"),
  alpha = 0.05
)

# 4. Apply apeglm log2 fold change shrinkage using the newly available coefficient
res_shrink <- lfcShrink(
  dds, 
  coef = "condition_pulldown_lab_vs_pulldown_unlab", 
  type = "apeglm"
)

# 4. View a quick summary of the results (How many genes are up/down regulated?)
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

# 1. Convert your shrunk results into a standard data frame
res_df <- as.data.frame(res_shrink) %>%
  rownames_to_column(var = "gene_id")

# ==============================================================================
# GENERATE TIER 1: Core High-Confidence Hits (The ~92 Genes)
# ==============================================================================
# Criteria: Statistically significant after multi-test correction, and enriched
tier1_core_hits <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  arrange(padj) # Sort by most significant

# Save Tier 1 to CSV
write.csv(tier1_core_hits, file = "pulldown_Tier1_core_hits_padj05.csv", row.names = FALSE)

# ==============================================================================
# GENERATE TIER 2: Broad Candidate Hits (The ~866 Genes)
# ==============================================================================
# Criteria: Significant by raw p-value, enriched, and excludes Tier 1 to avoid duplication
tier2_candidate_hits <- res_df %>%
  filter(pvalue < 0.05 & log2FoldChange > 0) %>%
  filter(!gene_id %in% tier1_core_hits$gene_id) %>% 
  arrange(pvalue) # Sort by raw p-value

# Save Tier 2 to CSV
write.csv(tier1_core_hits, file = "pulldown_Tier1_candidate_hits_pvalue05.csv", row.names = FALSE)
write.csv(tier2_candidate_hits, file = "pulldown_Tier2_candidate_hits_pvalue05.csv", row.names = FALSE)

# Print a quick confirmation to your console
message("Files successfully written!")
message("Tier 1 (Strict padj < 0.05): ", nrow(tier1_core_hits), " genes.")
message("Tier 2 (Raw pvalue < 0.05, excluding Tier 1): ", nrow(tier2_candidate_hits), " genes.") 

library(clusterProfiler)
library(org.Mm.eg.db)  # Loaded the mouse database
library(enrichplot)
library(dplyr)

# 1. Prepare your combined 866 gene list
target_genes <- c(tier1_core_hits$gene_id, tier2_candidate_hits$gene_id)

# 2. Run the Gene Ontology Enrichment for Mouse
ego_mouse <- enrichGO(gene          = target_genes,
                      OrgDb         = org.Mm.eg.db,  # Updated for mouse
                      keyType       = "SYMBOL",      # Assumes standard mouse symbols like 'Actb'
                      ont           = "BP",          # "BP" = Biological Process
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)

# 3. Save the mouse results table to a CSV file
go_mouse_results <- as.data.frame(ego_mouse)
write.csv(go_mouse_results, file = "mouse_pulldown_866_GO_results.csv", row.names = FALSE)

 # ==============================================================================
 # VISUALIZE THE MOUSE HITS
 # ==============================================================================

# # 1. Dotplot (Shows ratio of genes vs statistical significance)
 dotplot(ego_mouse, showCategory = 15) + ggtitle("Top 15 Mouse Biological Processes")

# # 2. Barplot (Shows absolute gene counts)
 barplot(ego_mouse, showCategory = 15) + ggtitle("Top 15 Mouse Biological Processes") 

# 1. Convert the GO object to a standard data frame if you haven't already
go_results <- as.data.frame(ego_mouse)

# 2. Filter for the specific ciliary terms
# This uses grepl to match "cilium assembly" or "cilium organization" case-insensitively
cilia_terms <- go_results %>%
  filter(grepl("cilium assembly|cilium organization", Description, ignore.case = TRUE))

# Print the terms found to make sure you have the right ones
print(cilia_terms[, c("ID", "Description", "Count")])

# 3. Extract and clean the gene list
# This grabs the 'geneID' column, splits the genes apart, and keeps only unique names
cilia_genes_list <- cilia_terms$geneID %>%
  paste(collapse = "/") %>%       # Combine if there are multiple matching rows
  strsplit(split = "/") %>%       # Split by the slash
  unlist() %>%                    # Flatten into a vector
  unique() %>%                    # Remove duplicates between the two terms
  sort()                          # Alphabetize

# Print the total count and the final list of genes to the console
message("\nFound ", length(cilia_genes_list), " unique genes matching cilia terms:")
print(cilia_genes_list)