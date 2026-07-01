
file_1 <- "o_outputs/processed_data/DE_results/DESeq2_upreg-Results_lab_pull_down_vs_unlab_pull_down.csv"
data_1 <- read.csv(file_1, stringsAsFactors = FALSE)

file_2 <- "o_outputs/processed_data/DE_results/DESeq2_upreg-Results_lab_pull_down_vs_lab_input.csv"
data_2 <- read.csv(file_2, stringsAsFactors = FALSE)

file_3 <- "o_outputs/processed_data/DE_results/DESeq2_upreg-Results_lab_pull_down_vs_unlab_input.csv"
data_3 <- read.csv(file_3, stringsAsFactors = FALSE)

file_4 <- "o_outputs/processed_data/DE_results/DESeq2_upreg-Results_unlab_pull_down_vs_lab_input.csv"
data_4 <- read.csv(file_4, stringsAsFactors = FALSE)

file_5 <- "o_outputs/processed_data/DE_results/DESeq2_upreg-Results_unlab_pull_down_vs_unlab_input.csv"
data_5 <- read.csv(file_5, stringsAsFactors = FALSE)

file_6 <- "o_outputs/processed_data/DE_results/DESeq2_upreg-Results_lab_input_vs_unlab_input.csv"
data_6 <- read.csv(file_6, stringsAsFactors = FALSE)

# for each entry in column 1 of file 1, check if it is present in column 1 of file 2, file 3, file 4, file 5, and file 6
for (i in seq_len(nrow(data_1))) {
    gene <- data_1[i, 1]
    data_1[i, "in_file_2"] <- gene %in% data_2[, 1]
    data_1[i, "in_file_3"] <- gene %in% data_3[, 1]
    data_1[i, "in_file_4"] <- gene %in% data_4[, 1]
    data_1[i, "in_file_5"] <- gene %in% data_5[, 1]
    data_1[i, "in_file_6"] <- gene %in% data_6[, 1]
}
