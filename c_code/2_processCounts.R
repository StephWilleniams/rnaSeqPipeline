
### Title: RNA-Seq Post-Counted Processing Workflow
### Authors: Emily Skates, Stephen Williams
## This script takes the output of the STAR + featureCounts, and performs a series of analyses

# Analysis includes:
# - PCA and PERMANOVA to assess sample clustering and group differences
# - Differential expression analysis using DESeq2
# - Co-expression network analysis using WGCNA
# - Gene Ontology (GO) enrichment analysis for up- and down-regulated genes

# Input libraries -- assumes they are installed, if not use install.packages("NAME") or BiocManager::install("NAME")
library(DESeq2)
library(ggplot2)
library(WGCNA)
library(pheatmap)
library(vegan)
library(pairwiseAdonis)
library(glue)
library(gtools)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)

# Set the comparison treatments for differential expression analysis
comparison_1 <- "lab_pull_down"
comparison_2 <- "unlab_pull_down"
copy_threshold <- 10 # Number of gene copys required to be considered expressed in a sample
condition_threshold <- 3 # Number of samples in a condition that must meet the copy threshold to be considered expressed in that condition
PCA_gene_count <- 500 # Number of top variable genes to include in the PCA analysis
network_size <- 100 # Number of top variable genes to include in the co-expression network analysis

# Run analysis options
PERFORM_DISPERSION_DIAGNOSTICS <- FALSE # Set to TRUE if you want to generate a dispersion diagnostic plot, may cause soft errors.
MAKE_PCA_PLOT <- FALSE                  # Set to TRUE if you want to generate a PCA plot.
DO_GROUP_COMPARISONS <- FALSE           # Set to TRUE if you want to perform PERMANOVA and Beta Diversity analysis.
DO_DE_ANALYSIS <- TRUE                 # Set to TRUE if you want to perform differential expression analysis.
PLOT_DE_ANALYSIS <- TRUE               # Set to TRUE if you want to generate a volcano plot for the differential expression analysis.
PLOT_MA_ANALYSIS <- TRUE               # Set to TRUE if you want to generate an MA plot for the differential expression analysis.
USE_LFC_SHRINK <- FALSE                 # Set to TRUE if you want to apply log fold change shrinkage to the DESeq2 results. This is optional and can help reduce noise in the results.
DO_GO_ANALYSIS <- TRUE                 # Set to TRUE if you want to perform Gene Ontology (GO) enrichment analysis on the up- and down-regulated genes.

# ====================================================================
# Phase 1a: Load Data & perform DESeq2 Setup
# ====================================================================

file_list <- Sys.glob("o_outputs/sample_*/gene_counts.txt")
file_list <- mixedsort(file_list)
temp_first <- read.table(file_list[1], header = TRUE, skip = 1, stringsAsFactors = FALSE)
count_matrix <- data.frame(Geneid = temp_first[, 1])
for (file in file_list) {
    sample_name <- basename(dirname(file))
    temp_data <- read.table(file, header = TRUE, skip = 1, stringsAsFactors = FALSE)
    count_matrix[[sample_name]] <- temp_data[, 7]
}
row.names(count_matrix) <- count_matrix$Geneid
count_matrix$Geneid <- NULL
sample_info <- data.frame(
    row.names = colnames(count_matrix),
    Condition1 = c(rep("unlab", 3), rep("lab", 3), rep("unlab", 3), rep("lab", 3)),
    Condition2 = c(rep("pull_down", 6), rep("input", 6))
)
sample_info$Group <- factor(paste0(sample_info$Condition1, "_", sample_info$Condition2))
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ Group)
keep <- rowSums(counts(dds) >= copy_threshold) >= condition_threshold
dds <- dds[keep, ]
vsd <- vst(dds, blind = FALSE) # Variance Stabilizing Transformation (VST)

# ====================================================================
# Phase 1b: Dispersion Diagnostics (OPTIONAL)
# ====================================================================

if (PERFORM_DISPERSION_DIAGNOSTICS) {

    dds <- dds |>
        estimateSizeFactors() |>
        estimateDispersions()
    disp_df <- data.frame(
        baseMean = mcols(dds)$baseMean,
        dispGeneEst = mcols(dds)$dispGeneEst,
        dispersion = mcols(dds)$dispersion,
        dispFit = mcols(dds)$dispFit
    )
    disp_df <- disp_df[disp_df$baseMean > 0, ]

    dispersion_plot <- ggplot(disp_df, aes(x = baseMean)) +
        geom_point(aes(y = dispGeneEst), color = "black", alpha = 0.3, size = 1) + # Plot raw gene estimates (equivalent to black dots)
        geom_point(aes(y = dispersion), color = "dodgerblue", alpha = 0.3, size = 1, shape = 1) + # Plot final shrink/MAP estimates (equivalent to blue circles overlay)
        geom_line(aes(y = dispFit), color = "red", size = 1) + # Plot the fitted dispersion trend line (equivalent to the red curve)
        scale_x_log10() +
        scale_y_log10(limits = c(1e-4, 10)) +
        theme_minimal() +
        labs(title = "DESeq2 Dispersion Estimates Diagnostics",
             x = "Mean of Normalized Counts",
             y = "Dispersion")

    ggsave(filename = glue("f_figures/Dispersion_Estimates_{comparison_1}_vs_{comparison_2}.png"),
           plot = dispersion_plot,
           width = 8,
           height = 6,
           dpi = 300)
}

# ====================================================================
# Phase 2a: PCA (Checking how your groups cluster)
# ====================================================================

if (MAKE_PCA_PLOT) {

    pca_data <- plotPCA(vsd, intgroup = c("Group"), ntop = PCA_gene_count, returnData = TRUE)

    pca_figure <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, shape = group)) +
        geom_point(size = 2.5) + 
        theme_minimal() +

        scale_color_manual(
            name = "Experimental Condition",
            values = c(
                "lab_input"        = "#D55E00", 
                "lab_pull_down"    = "#009E73", 
                "unlab_input"      = "#56B4E9", 
                "unlab_pull_down"  = "#CC79A7"  
            ),
            labels = c(
                "lab_input"        = "Labeled Input",
                "lab_pull_down"    = "Labeled Pull-Down",
                "unlab_input"      = "Unlabeled Input",
                "unlab_pull_down"  = "Unlabeled Pull-Down"
            )
        ) +
        scale_shape_manual(
            name = "Experimental Condition",
            values = c(
                "lab_input"        = 16,
                "lab_pull_down"    = 17,
                "unlab_input"      = 15,
                "unlab_pull_down"  = 18  
            ),
            labels = c(
                "lab_input"        = "Labeled Input",
                "lab_pull_down"    = "Labeled Pull-Down",
                "unlab_input"      = "Unlabeled Input",
                "unlab_pull_down"  = "Unlabeled Pull-Down"
            )
        ) +
        labs(
            x = paste0("PC1: ", round(attr(pca_data, "percentVar")[1] * 100, 1), "% variance"),
            y = paste0("PC2: ", round(attr(pca_data, "percentVar")[2] * 100, 1), "% variance")
        ) +
        theme(
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 14),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
        )

    ggsave(filename = glue("f_figures/PCA_Plot.png"),
           plot = pca_figure,
           width = 8,
           height = 4,
           dpi = 300)
}

# ====================================================================
# Phase 2b: PERMANOVA + Beta Diversity (Testing if your groups are significantly different)
# ====================================================================

if (DO_GROUP_COMPARISONS) {

    # Format the data for analysis
    mat_vst <- t(assay(vsd))
    dist_matrix <- vegdist(mat_vst, method = "euclidean")

    # PERMANOVA between 4 catergories (unlab_pull_down, lab_pull_down, unlab_input, lab_input)
    permanova_results <- adonis2(dist_matrix ~ Group, data = sample_info, permutations = 999)
    print("--------------------------------")
    print(permanova_results)
    print("--------------------------------")

    # PERMANOVA between 2 catergories (pull_down, input)
    permanova_results <- adonis2(dist_matrix ~ Condition2, data = sample_info, permutations = 999)
    print(permanova_results)
    print("--------------------------------")

    # # Pairwise PERMANOVA
    # pairwise_results <- pairwise.adonis(dist_matrix,
    #                                     factors = sample_info$Group,
    #                                     sim.method = "euclidean",
    #                                     p.adjust.m = "bonferroni") # Bonferroni is best for small n=12
    # print(pairwise_results)

    # Beta Dispersion Test
    dispersion_test <- betadisper(dist_matrix, group = sample_info$Group)
    anova_results <- anova(dispersion_test)
    print(anova_results)
    print("--------------------------------")
}

# ====================================================================
# Phase 3a: Differential Expression Analysis
# ====================================================================

if (DO_DE_ANALYSIS) {

    # dds <- estimateSizeFactors(dds) # DEFAULT SETTING
    # sizeFactors(dds) <- rep(1, ncol(dds)) # IGNORE DATASET SIZE (NOT RECOMMENDED)
    dds <- DESeq(dds)
    res_unshrunk <- results(dds, contrast = c("Group", comparison_1, comparison_2))

    if (USE_LFC_SHRINK) {
        res <- lfcShrink(dds,
                         contrast = c("Group", comparison_1, comparison_2),
                         res = res_unshrunk,
                         type = "ashr") # type="ashr" allows you to use the exact same 'contrast' argument
    } else {
        res <- res_unshrunk
    }

    upregs <- res[which(res$log2FoldChange > 0.5 & res$padj < 0.05), ]
    downregs <- res[which(res$log2FoldChange < -0.5 & res$padj < 0.05), ]

    print(paste("Number of Up-regulated Genes:", nrow(upregs)))
    print(paste("Number of Down-regulated Genes:", nrow(downregs)))
    write.csv(as.data.frame(upregs), file = paste0("o_outputs/processed_data/DE_results/DESeq2_upreg-Results_", comparison_1, "_vs_", comparison_2, ".csv"))
    write.csv(as.data.frame(downregs), file = paste0("o_outputs/processed_data/DE_results/DESeq2_downreg-Results_", comparison_1, "_vs_", comparison_2, ".csv"))
}

# ====================================================================
# Phase 3b: Volcano Plot (Visualizing Fold Change vs. Significance)
# ====================================================================

if (DO_DE_ANALYSIS && PLOT_DE_ANALYSIS) {

    # 1. Convert the DESeq2 results object into a standard R data frame
    res_df <- as.data.frame(res)

    # 2. Remove any rows with NA values in the adjusted p-value column
    res_df <- res_df[!is.na(res_df$padj), ] # (DESeq2 sets NA for genes with extremely low counts to save processing time)

    # 3. Create a new column to classify genes for coloring
    res_df$Significance <- "Not Significant"
    res_df$Significance[res_df$log2FoldChange > 0.5 & res_df$padj < 0.05] <- "Up-regulated"
    res_df$Significance[res_df$log2FoldChange < -0.5 & res_df$padj < 0.05] <- "Down-regulated"

    # 4. Build the plot using ggplot2
    volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
        geom_point(alpha = 0.6, size = 1.5) + # alpha adds slight transparency to see overlapping points
        scale_color_manual(values = c("Up-regulated" = "red",
                                      "Down-regulated" = "blue",
                                      "Not Significant" = "grey")) +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        theme_minimal() +
        labs(x = "Log2 Fold Change",
             y = "-Log10(Adjusted P-value)") +
        theme(legend.position = "right")

    ggsave(filename = paste0("f_figures/DE/Volcano_Plot_", comparison_1, "_vs_", comparison_2, ".png"),
           plot = volcano_plot,
           width = 8,
           height = 4,
           dpi = 300)
} else if (!DO_DE_ANALYSIS && PLOT_DE_ANALYSIS) {
    print("Differential expression analysis was not performed, no volcano plot can be generated.")
}

# ====================================================================
# Phase 3c: MA-Plot (Visualizing Mean Expression vs. Fold Change)
# ====================================================================

if (DO_DE_ANALYSIS && PLOT_MA_ANALYSIS) {
    png(glue("f_figures/MA/MA_Plot_{comparison_1}_vs_{comparison_2}.png"), width = 2400, height = 1800, res = 300)
    plotMA(res, main = glue("MA Plot: {comparison_1} vs {comparison_2}"), ylim = c(-4, 4))
    dev.off()

} else if (!DO_DE_ANALYSIS && PLOT_DE_ANALYSIS) {
    print("Differential expression analysis was not performed, no MA plot can be generated.")
}

# ====================================================================
# Phase 4: Network Analysis (Gene Correlation)
# ====================================================================

if (DO_DE_ANALYSIS) {

    norm_counts <- assay(vsd)
    dat_expr <- t(norm_counts)
    gene_correlations <- cor(dat_expr, method = "pearson")
    top_var_genes <- head(order(rowVars(norm_counts), decreasing = TRUE), network_size)

    pheatmap(norm_counts[top_var_genes, ],
             annotation_col = sample_info,
             show_rownames = FALSE,
             main = paste0("Co-expression of Top ", network_size, " Variable Genes"),
             filename = glue("f_figures/Network-coexpr/Top_", network_size, "_Coexpression_Heatmap_{comparison_1}_vs_{comparison_2}.png"),
             width = 8,
             height = 6)
}

# ====================================================================
# Phase 5a: Gene onology (GO) Enrichment Analysis
# ====================================================================

if (DO_DE_ANALYSIS && DO_GO_ANALYSIS) {

    up_genes <- rownames(upregs)
    down_genes <- rownames(downregs)
    universe_genes <- rownames(dds)
    org_db <- org.Mm.eg.db
    id_type <- "SYMBOL"

    go_up <- enrichGO(gene          = up_genes,
                      #universe      = universe_genes,
                      OrgDb         = org_db,
                      keyType       = id_type,
                      ont           = "ALL", # "BP" = Biological Process, "CC" = Cellular Component, "MF" = Molecular Function, "ALL" = All three
                      pAdjustMethod = "BH", # Benjamini-Hochberg false discovery rate
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE) # Converts IDs back to user-friendly gene symbols

    go_down <- enrichGO(gene        = down_genes,
                        # universe    = universe_genes,
                        OrgDb       = org_db,
                        keyType     = id_type,
                        ont         = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.2,
                        readable      = TRUE)

    write.csv(as.data.frame(go_up), file = glue("o_outputs/processed_data/GO_results/GO_Enrichment_UP_{comparison_1}_vs_{comparison_2}.csv"))
    write.csv(as.data.frame(go_down), file = glue("o_outputs/processed_data/GO_results/GO_Enrichment_DOWN_{comparison_1}_vs_{comparison_2}.csv"))
} else if (!DO_DE_ANALYSIS && DO_GO_ANALYSIS) {
    print("Differential expression analysis was not performed, no GO analysis can be generated.")
}

# ====================================================================
# Phase 5b: Gene onology (GO) Enrichment Plot
# ====================================================================

if (DO_DE_ANALYSIS && DO_GO_ANALYSIS && PLOT_DE_ANALYSIS) {

    if (!is.null(go_up) && nrow(go_up) > 0) {

        dot_up <- dotplot(go_up, showCategory = 20) +
            ggtitle(glue("Top Biological Processes: Up-regulated in {comparison_1}"))

        ggsave(filename = glue("f_figures/GO/GO_Dotplot_UP_{comparison_1}_vs_{comparison_2}.png"),
               plot = dot_up, width = 9, height = 7, dpi = 300)

        # # Barplot: standard enrichment visual
        # bar_up <- barplot(go_up) +
        #     ggtitle(paste0("Top Biological Processes: Up-regulated in ", comparison_1))

        # ggsave(filename = glue("f_figures/GO_Barplot_UP_{comparison_1}_vs_{comparison_2}.png"),
        #        plot = bar_up, width = 9, height = 7, dpi = 300)
    }

    if (!is.null(go_down) && nrow(go_down) > 0) {

        dot_down <- dotplot(go_down, showCategory = 20) +
            ggtitle(glue("Top Biological Processes: Down-regulated in {comparison_1}"))

        ggsave(filename = glue("f_figures/GO/GO_Dotplot_DOWN_{comparison_1}_vs_{comparison_2}.png"),
               plot = dot_down, width = 9, height = 7, dpi = 300)

        # bar_down <- barplot(go_down, showCategory = 20) +
        #     ggtitle(glue("Top Biological Processes: Down-regulated in {comparison_1}"))

        # ggsave(filename = glue("f_figures/GO_Barplot_DOWN_{comparison_1}_vs_{comparison_2}.png"),
        #        plot = bar_down, width = 9, height = 7, dpi = 300)
    }
} else if ((!DO_DE_ANALYSIS && DO_GO_ANALYSIS && PLOT_DE_ANALYSIS) || (DO_DE_ANALYSIS && !DO_GO_ANALYSIS && PLOT_DE_ANALYSIS)) {
    print("Required analysis was not performed, no GO plots can be generated.")
}
