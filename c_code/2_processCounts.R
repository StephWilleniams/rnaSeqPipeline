
# Install required packages if missing
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "WGCNA", "pheatmap"))

# Input libraries -- assumes they are installed, if not use install.packages("NAME") or BiocManager::install("NAME")
library(DESeq2)
library(ggplot2)
library(WGCNA)
library(pheatmap)
library(vegan)
library(pairwiseAdonis)
library(glue)

comparison_1 <- "lab_input"
comparison_2 <- "unlab_input"

# ====================================================================
# Phase 1a: Preparing Your Data & DESeq2 Setup
# ====================================================================

# 1. Find all 12 'gene_counts.txt' files dynamically
# Make sure this matches the OUT_DIR structure from your bash script
file_list <- Sys.glob("o_outputs/sample_*/gene_counts.txt")
# print(paste("Found", length(file_list), "files."))

# Read the first file to extract the Gene IDs
temp_first <- read.table(file_list[1], header = TRUE, skip = 1, stringsAsFactors = FALSE)
count_matrix <- data.frame(Geneid = temp_first[, 1])

# Loop through the matched files and grab ONLY the counts (Column 7)
for (file in file_list) {
    sample_name <- basename(dirname(file)) # Extracts "sample_X"
    temp_data <- read.table(file, header = TRUE, skip = 1, stringsAsFactors = FALSE)
    count_matrix[[sample_name]] <- temp_data[, 7]
}

# Clean up the matrix formatting
row.names(count_matrix) <- count_matrix$Geneid
count_matrix$Geneid <- NULL

# 2. Create your metadata table (colData)
sample_info <- data.frame(
    row.names = colnames(count_matrix),
    Condition1 = c(rep("unlab", 3), rep("lab", 3), rep("unlab", 3), rep("lab", 3)),
    Condition2 = c(rep("pull_down", 6), rep("input", 6))
)

sample_info$Group <- factor(paste0(sample_info$Condition1, "_", sample_info$Condition2))

# 3. Build the DESeq2 Object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ Group)

copy_threshold <- 10
condition_threshold <- 3

keep <- rowSums(counts(dds) >= copy_threshold) >= condition_threshold
dds <- dds[keep, ]

# ====================================================================
# Phase 1b: Dispersion Diagnostics (ggplot version)
# ====================================================================

# 1. To get dispersion estimates before the full DESeq run, we run these two steps:
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# 2. Extract dispersion metrics into a dataframe for ggplot
disp_df <- data.frame(
    baseMean = mcols(dds)$baseMean,
    dispGeneEst = mcols(dds)$dispGeneEst,
    dispersion = mcols(dds)$dispersion,
    dispFit = mcols(dds)$dispFit
)

# Filter out rows with 0 baseMean to prevent log10 errors
disp_df <- disp_df[disp_df$baseMean > 0, ]

# 3. Create the plot
dispersion_plot <- ggplot(disp_df, aes(x = baseMean)) +
    # Plot raw gene estimates (equivalent to black dots)
    geom_point(aes(y = dispGeneEst), color = "black", alpha = 0.3, size = 1) +
    # Plot final shrink/MAP estimates (equivalent to blue circles overlay)
    geom_point(aes(y = dispersion), color = "dodgerblue", alpha = 0.3, size = 1, shape = 1) +
    # Plot the fitted dispersion trend line (equivalent to the red curve)
    geom_line(aes(y = dispFit), color = "red", size = 1) +
    # Apply standard RNA-seq log scales to both axes
    scale_x_log10() +
    scale_y_log10(limits = c(1e-4, 10)) +
    theme_minimal() +
    labs(title = "DESeq2 Dispersion Estimates Diagnostics",
         x = "Mean of Normalized Counts",
         y = "Dispersion")

# Print to console/RStudio viewer
print(dispersion_plot)

# 4. Save the figure natively via ggsave
ggsave(filename = glue("f_figures/Dispersion_Estimates_{comparison_1}_vs_{comparison_2}.png"),
       plot = dispersion_plot,
       width = 8,
       height = 6,
       dpi = 300)

# ====================================================================
# Phase 2a: PCA (Checking how your ND groups cluster)
# ====================================================================

# Variance Stabilizing Transformation (VST)
vsd <- vst(dds, blind = FALSE)

# Assign the plot to a variable named 'pca_figure'
pca_figure <- plotPCA(vsd, intgroup = c("Group")) +
    # stat_ellipse(aes(color = Group), type = "t") +
    theme_minimal() +
    ggtitle("PCA of 12 Samples: treatment Groups vs Cell Parts")

# SAVE THE FIGURE
# Tell ggsave to save the 'pca_figure' variable we just created
ggsave(filename = glue("f_figures/PCA_Plot.png"),
       plot = pca_figure,
       width = 8,
       height = 6,
       dpi = 300)

# ====================================================================
# Phase 2b: PERMANOVA + Beta Diversity (Testing if your groups are significantly different)
# ====================================================================

# 1. Extract your normalized counts matrix from the VST object
# We transpose it (t) because vegan expects samples as rows and genes as columns
mat_vst <- t(assay(vsd))

# 2. Calculate a Euclidean Distance Matrix among all 12 samples
dist_matrix <- vegdist(mat_vst, method = "euclidean")

# 3. Run the PERMANOVA
# This asks: "Does our 'Group' column significantly separate the distances?"
permanova_results <- adonis2(dist_matrix ~ Group, data = sample_info)

# Print the table to see your P-value and R2 effect size!
print(permanova_results)

# 2. Perform pairwise comparisons using the Euclidean distance matrix
# This runs the adonis test for every combination of your 'Group' column
pairwise_results <- pairwise.adonis(dist_matrix,
                                    factors = sample_info$Group,
                                    sim.method = "euclidean",
                                    p.adjust.m = "bonferroni") # Bonferroni is best for small n=12

# 3. View the table of significant differences
print(pairwise_results)

# Calculate multivariate homogeneity of groups dispersals
dispersion_test <- betadisper(dist_matrix, group = sample_info$Group)

# Run a permutation test to see if group spreads are statistically identical
anova_results <- anova(dispersion_test)

# Print the table to see your P-value and R2 effect size!
print(anova_results)

# ====================================================================
# Phase 3a: Differential Expression Analysis
# ====================================================================

# THINGS TO TRY


# 1. Run the core DESeq2 algorithm (This handles your library depth normalization!)

# sizeFactors(dds) <- rep(1, ncol(dds))
dds <- DESeq(dds)

# 2. Extract the un-shrunk results (Standard)
res_unshrunk <- results(dds, contrast = c("Group", comparison_1, comparison_2))

# 3. Apply lfcShrink to clean up the noise
# type="ashr" allows you to use the exact same 'contrast' argument
res <- lfcShrink(dds,
                 contrast = c("Group", comparison_1, comparison_2),
                 res = res_unshrunk,
                 type = "ashr")

# List the significant genes (p-value < 0.05) and sort by log2 fold change
upregs <- res_unshrunk[which(res_unshrunk$log2FoldChange > 0.5 & res_unshrunk$pvalue < 0.05), ]
downregs <- res_unshrunk[which(res_unshrunk$log2FoldChange < -0.5 & res_unshrunk$pvalue < 0.05), ]

print(paste("Number of Up-regulated Genes:", nrow(upregs)))
print(paste("Number of Down-regulated Genes:", nrow(downregs)))

# Save the results to a CSV file for later use
write.csv(as.data.frame(upregs), file = glue("o_outputs/processed_data/DESeq2_upreg-Results_{comparison_1}_vs_{comparison_2}.csv"))
write.csv(as.data.frame(downregs), file = glue("o_outputs/processed_data/DESeq2_downreg-Results_{comparison_1}_vs_{comparison_2}.csv"))

# ====================================================================
# Phase 3b: Volcano Plot (Visualizing Fold Change vs. Significance)
# ====================================================================

# 1. Convert the DESeq2 results object into a standard R data frame
res_df <- as.data.frame(res_unshrunk)

# 2. Remove any rows with NA values in the adjusted p-value column
# (DESeq2 sets NA for genes with extremely low counts to save processing time)
res_df <- res_df[!is.na(res_df$pvalue), ]

# 3. Create a new column to classify genes for coloring
# We will flag genes as "Significant" if their adjusted p-value is < 0.05
# AND their Log2 Fold Change is greater than 0.5 or less than -0.5 (a 1.5-fold change).
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$log2FoldChange > 0.5 & res_df$pvalue < 0.05] <- "Up-regulated"
res_df$Significance[res_df$log2FoldChange < -0.5 & res_df$pvalue < 0.05] <- "Down-regulated"

# 4. Build the plot using ggplot2
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) + # alpha adds slight transparency to see overlapping points

    # Set custom colors for the three groups
    scale_color_manual(values = c("Up-regulated" = "red",
                                  "Down-regulated" = "blue",
                                  "Not Significant" = "grey")) +

    # Add dashed lines to show our thresholds
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +

    # Clean up the theme and labels
    theme_minimal() +
    labs(title = "Volcano Plot: Labeled {comparison_1} vs Labeled {comparison_2}",
         x = "Log2 Fold Change",
         y = "-Log10(Adjusted P-value)") +
    theme(legend.position = "right")

# Display the plot in your R session
# print(volcano_plot)

# 5. Save the figure (just like we did for the PCA!)
ggsave(filename = glue("f_figures/Volcano_Plot_{comparison_1}_vs_{comparison_2}.png"),
       plot = volcano_plot,
       width = 8,
       height = 6,
       dpi = 300)

# ====================================================================
# Phase 3c: MA-Plot (Visualizing Mean Expression vs. Fold Change)
# ====================================================================

# 1. Display the plot in your RStudio viewer
# ylim restricts the y-axis so massive outliers don't squash the rest of the data
plotMA(res, main = glue("MA Plot: {comparison_1} vs {comparison_2}"), ylim = c(-4, 4))

# 2. Save the plot to your figures folder
png(glue("f_figures/MA_Plot_{comparison_1}_vs_{comparison_2}.png"), width = 2400, height = 1800, res = 300)

# Redraw the plot inside the PNG file
plotMA(res, main = glue("MA Plot: {comparison_1} vs {comparison_2}"), ylim = c(-4, 4))

# Close and save the file
dev.off()

# ====================================================================
# Phase 4: Network Analysis (Gene Correlation)
# ====================================================================

# Extract normalized counts for network analysis
norm_counts <- assay(vsd)
dat_expr <- t(norm_counts)

# Calculate the correlation matrix (base for WGCNA)
gene_correlations <- cor(dat_expr, method = "pearson")

# Top 1000 most variable genes for the heatmap
top_var_genes <- head(order(rowVars(norm_counts), decreasing = TRUE), 1000)

# Plot AND Save Co-expression Heatmap
# CORRECTED: Changed 'topVarGenes' to 'top_var_genes' to match the variable name
pheatmap(norm_counts[top_var_genes, ],
         annotation_col = sample_info,
         show_rownames = FALSE,
         main = "Co-expression of Top 1000 Variable Genes",
         filename = glue("f_figures/Top_1000_Coexpression_Heatmap_{comparison_1}_vs_{comparison_2}.png"),
         width = 8,
         height = 6)