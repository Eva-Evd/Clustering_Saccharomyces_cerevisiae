rm(list = ls())
################################################################################
################################################################################
# import data
old_data <- read.csv("generations.csv")

# View the first few rows
head(old_data)
################################################################################
################################################################################
sum(is.na(old_data))

library(tidyverse)

# 1. Load the 2026 significant genes list
sig_genes_2026 <- read.csv("significant_aging_genes_2026.csv")

# 2. Clean up identifiers, drop NAs, AND filter by the 2026 gene list
cleaned_aging_data <- old_data %>%
  filter(gene != "" & gene != "-") %>%
  # Extract primary systematic/standard name if separated by pipes
  mutate(gene = str_split_i(gene, "\\|", 1)) %>%
  # CRUCIAL STEP: Keep ONLY genes that were found significant in the 2026 study
  filter(gene %in% sig_genes_2026$gene) %>%
  # If a gene appears twice, average its values
  group_by(gene) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  ungroup() %>%
  # Drop any rows that still contain NA values in any column
  drop_na()

cat("Number of validated cross-study genes remaining for clustering:", nrow(cleaned_aging_data), "\n")



# 1. Build the mathematical matrix (drop the first text column 'gene')
expr_matrix <- as.matrix(cleaned_aging_data[, -1])
rownames(expr_matrix) <- cleaned_aging_data$gene

# 2. Z-score Scale the rows so genes can be compared by trend rather than raw value
scaled_aging <- t(scale(t(expr_matrix)))

################################################################################
################################################################################
# Calculate Pearson correlation matrix across genes
gene_cor <- cor(t(scaled_aging), method = "pearson")

# Convert correlation to a distance metric (bounded between 0 and 2)
gene_dist <- as.dist(1 - gene_cor)

# Cluster 
hc_results <- hclust(gene_dist, method = "ward.D2")



# Based on the paper's findings, let's cut into 6 distinct aging profiles
# (e.g., Early stress response, late down-regulation, stationary phase, etc.)
clusters <- cutree(hc_results, k = 6)
# Extract cluster sizes
cluster_sizes <- table(clusters)

# View the result
print(cluster_sizes)

################################################################################
################################################################################

# Bind cluster info to data and pivot to long-format for ggplot
aging_profiles <- as.data.frame(scaled_aging) %>%
  rownames_to_column(var = "gene") %>%
  mutate(Cluster = as.factor(clusters)) %>%
  pivot_longer(cols = starts_with("gen_"), 
               names_to = "Generation", 
               values_to = "Relative_Expression") %>%
  # Order the X-axis chronologically
  mutate(Generation = factor(Generation, levels = c("gen_1", "gen_8", "gen_12", "gen_19")))

# Plot the expression tracks
ggplot(aging_profiles, aes(x = Generation, y = Relative_Expression, group = gene)) +
  # Light blue background lines for individual yeast genes
  geom_line(alpha = 0.15, color = "skyblue") + 
  # Bold red line tracking the overall trend of that aging cluster
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "firebrick", linewidth = 1.3) +
  facet_wrap(~ Cluster, labeller = label_both) +
  theme_bw() +
  labs(title = "S. cerevisiae Replicative Aging Expression Profiles (GSE10018)",
       subtitle = "Pearson Correlation Distance Clustering",
       y = "Relative Expression (Z-score)",
       x = "Replicative Age (Generations)")




#########################################################################################################
#########################################################################################################
################################################################################
#  Standalone Colored Dendrogram
################################################################################
library(dendextend)

# 1. Convert the hclust object into a formal dendrogram object
dend <- as.dendrogram(hc_results)

# 2. Color the branches dynamically based on your 3 clusters
# (Using clean, high-contrast colors)
dend <- color_branches(dend, k = 3, col = c("firebrick3", "navy", "darkorange2", "forestgreen"))

# 3. Clean up the label formatting (hiding individual names since 164 is too crowded)
dend <- set(dend, "labels_cex", 0.0001)

# 4. Plot the final high-resolution dendrogram
par(mar = c(3, 5, 4, 2)) # Adjust plot margins
plot(dend, 
     main = "Hierarchical Clustering Tree of 448 Core Aging Genes",
     ylab = "Distance (1 - Pearson Correlation)",
     sub = "Cut into 4 distinct expression trajectories",
     leaflab = "none") # Completely hides individual gene text overlap

# 5. Add a horizontal line demonstrating exactly where cutree() sliced the branches
rect.dendrogram(dend, k = 3, border = "gray60", lty = 2, lwd = 1.5)


################################################################################
#  Extract and Export Gene Lists per Cluster (PathBIX  Optimized)
################################################################################
cat("\n--- Exporting PathBIX-Ready Gene Lists to Working Directory ---\n")

for(i in 1:3) {
  # 1. Isolate the names of the genes falling into cluster i
  cluster_genes <- names(clusters[clusters == i])
  
  # 2. PathBIX Clean-Up: Remove any background artifact characters or punctuation
  # Ensures names look like standard systematic (e.g., YAL018C) or standard common names (e.g., ETS1-1)
  cluster_genes <- trimws(cluster_genes)               # Strip hidden spaces
  cluster_genes <- cluster_genes[cluster_genes != ""]  # Drop empty rows
  cluster_genes <- unique(cluster_genes)               # Ensure absolute uniqueness per file
  
  # 3. Create a clean data frame with NO column header line
  # PathBIX background parsers prefer a raw text/csv newline list without a "gene" header row
  export_df <- data.frame(gene = cluster_genes)
  
  # 4. Define a unique filename
  filename <- paste0("pathbix_cluster_", i, "_genes.txt")
  
  # 5. Write out as a clean, headerless plain text list (one gene per line)
  # col.names = FALSE removes the text "gene" header so PathBIX reads it instantly.
  # quote = FALSE removes quotation marks so the raw string is cleanly parsed.
  write.table(export_df, 
              file = filename, 
              row.names = FALSE, 
              col.names = FALSE, 
              quote = FALSE, 
              sep = "\n")
  
  cat("Saved PathBIX Cluster", i, "(", length(cluster_genes), " genes) -> ", filename, "\n")
}