rm(list = ls())
################################################################################
# 1. Import and Inspect Data
################################################################################
counts <- read.csv("wt_data.csv")

library(edgeR)
library(pheatmap)
library(DESeq2)
library(dplyr)
library(ggplot2)

dim(counts)
head(counts)
sum(is.na(counts))

################################################################################
# 2. Correct Pre-Processing: Establish Row Names But Keep Original Object
################################################################################
# Step A: Move gene_id to row names on a numeric-only matrix so the math works
counts_numeric <- counts[, -1] 
rownames(counts_numeric) <- counts$gene_id

# Step B: Round decimals to integers so DESeq2/edgeR won't crash
counts_numeric <- round(counts_numeric)
head(counts_numeric) 

# 3. Plot Library Sizes
lib_sizes <- colSums(counts_numeric)
labels <- names(lib_sizes)
colors <- c(rep("skyblue", 3), rep("darkorange", 3))

par(mar = c(7, 5, 4, 2)) 
barplot(lib_sizes,
        col = colors,
        names.arg = labels,
        las = 2,           
        cex.names = 0.8,    
        main = "Library Sizes per Sample",
        ylab = "Total Reads / Counts")

################################################################################
# 4. Filter Out Low Counts (Safely Keeping Gene IDs)
################################################################################
sample_names <- colnames(counts)[-1]
time_points <- sapply(strsplit(sample_names, "_"), `[`, 2)
time_points <- gsub("log", "0", time_points)
groups <- factor(time_points, levels = c("0", "24"))

# Run filterByExpr on the numeric-only matrix
keep <- filterByExpr(counts_numeric, group = groups)

# Subset your original data frame to KEEP the gene_id column intact
counts_filtered <- counts[keep, ]

nrow(counts_filtered)
head(counts_filtered)

################################################################################
# 5. Format Matrix and Run DESeq2
################################################################################
# Drop the gene_id column for the mathematical matrix, but assign row names
counts_matrix <- as.matrix(counts_filtered[, -1])  
rownames(counts_matrix) <- counts_filtered$gene_id 

# Ensure whole integers
counts_matrix <- round(counts_matrix)

# Create metadata
metadata <- data.frame(
  Time = factor(c("0", "0", "0", "24", "24", "24"), levels = c("0", "24")),
  row.names = colnames(counts_matrix)
)

# Create and run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = metadata,
  design = ~ Time
)

dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)
mat <- assay(vsd)

################################################################################
# 6. Quality Control: PCA Plot
################################################################################
resultsNames(dds)
dim(mat)
summary(as.vector(mat))

pca_data <- plotPCA(vsd, intgroup = "Time", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = Time)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  scale_color_manual(values = c("0" = "blue", "24" = "red")) + 
  theme_minimal(base_size = 14) +
  labs(title = "PCA Plot", color = "Time Point (Hours)")

################################################################################
# 7. Gene Selection & Extraction (Fixing the Missing Lists Error)
################################################################################
res <- results(dds, name = "Time_24_vs_0")
res_df <- as.data.frame(res)

alpha <- 0.05
lfc_threshold <- 1

# Filter for all significant genes
sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < alpha & abs(res_df$log2FoldChange) > lfc_threshold, ]
cat("Total significant genes:", nrow(sig_genes), "\n")

# Extract the universal list of significant gene IDs
significant_gene_ids <- rownames(sig_genes)




# Create a data frame containing just your significant gene IDs
gene_export <- data.frame(gene = significant_gene_ids)

# Export to a CSV file (no row names, so it stays clean)
write.csv(gene_export, file = "significant_aging_genes_2026.csv", row.names = FALSE)

cat("Successfully exported", nrow(gene_export), "genes to CSV!\n")

#################################################################################













# 1. Load the required annotation library
# (If not installed, run: BiocManager::install("org.Sc.sgd.db"))
library(org.Sc.sgd.db)

# 2. Read your current file containing ORF names
sig_genes_2026 <- read.csv("significant_aging_genes_2026.csv")

# 3. Extract the ORF list (assuming the column name is "gene" or "ORF")
# Adjust the '$gene' part if your column has a different header name
orf_keys <- sig_genes_2026$gene

# 4. Map ORF names to Standard Common Gene Names
gene_symbols <- mapIds(
    org.Sc.sgd.db,
    keys = orf_keys,
    column = "GENENAME",   # This extracts common symbols like 'ETS1-1'
    keytype = "ORF",       # This specifies your input type is systematic 'YAL018C'
    multiVals = "first"    # If an ORF matches multiple symbols, take the first one
)

# 5. Create a new data frame using the common gene symbols
updated_gene_export <- data.frame(gene = gene_symbols)

# 6. Optional Clean-up: If an ORF doesn't have a common symbol, R returns NA.
# If you want to keep the original ORF name instead of a blank/NA, run this line:
updated_gene_export$gene[is.na(updated_gene_export$gene)] <- orf_keys[is.na(updated_gene_export$gene)]

# If you prefer to completely remove genes that don't have common names, run this instead:
# updated_gene_export <- updated_gene_export[!is.na(updated_gene_export$gene), , drop = FALSE]

# 7. Overwrite the CSV file with the updated common names
write.csv(updated_gene_export, file = "significant_aging_genes_2026.csv", row.names = FALSE)

cat("Successfully converted ORF IDs to Common Gene Names and saved the file!\n")