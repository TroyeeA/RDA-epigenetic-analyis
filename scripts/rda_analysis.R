################################################################################
# Redundancy Analysis (RDA) for Genotype-Environment Association
################################################################################
# 
# Description: RDA analysis to determine how groups of loci covary in response 
#              to multivariate environments. This script performs RDA on genetic
#              (SNP) and epigenetic (methylation) data.
#
# Author: A. N. Troyee
# Date: 2024
# License: MIT
#
# References:
#   - Legendre P, Legendre L (2012) Numerical Ecology, 3rd edition. Elsevier.
#   - Borcard D, Gillet F, Legendre P (2011) Numerical Ecology with R. Springer.
#   - Oksanen J et al. (2016) vegan: Community Ecology Package.
#   - https://popgen.nescent.org/2018-03-27_RDA_GEA.html
#   - https://www.davidzeleny.net/anadat-r/doku.php/en:rda_cca
#
# Important Note:
#   RDA is a linear model and assumes linear dependence between response 
#   variables (genotypes) and explanatory variables (environmental predictors).
################################################################################

# Load required packages -------------------------------------------------------
suppressPackageStartupMessages({
  library("data.table")
  library("ggplot2")
  library("tidyverse")
  library("gridExtra")
  library("hierfstat")  # Compute F statistics
  library("impute")
  library("psych")      # Investigate correlations among predictors
  library("vegan")      # Run RDA
  library("adegenet")   # Get allele counts
})

# Set working directory (modify as needed) -------------------------------------
# setwd("path/to/your/project")

# Define input file paths ------------------------------------------------------
# Modify these paths to match your file locations
input_dir <- "data"
output_dir <- "results"
script_dir <- "scripts"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Input files
sample_table_file <- file.path(input_dir, "sample_design.txt")
methylation_file <- file.path(input_dir, "methylation_data.bed")
snp_file <- file.path(input_dir, "snp_data.txt")

# Optional: Load custom functions if available
# source(file.path(script_dir, "custom_functions.R"))

################################################################################
# PART 1: LOAD AND PREPARE DATA
################################################################################

# Load sample metadata ---------------------------------------------------------
cat("Loading sample metadata...\n")
sample_table <- read.table(sample_table_file, 
                           sep = '\t', 
                           header = TRUE, 
                           stringsAsFactors = FALSE, 
                           row.names = 1)

sample_ids <- rownames(sample_table)

# Load SNP data ----------------------------------------------------------------
cat("Loading SNP data...\n")
# Assuming SNP file is tab-delimited with samples as columns
snp_data <- read.table(snp_file, 
                       sep = '\t', 
                       header = TRUE, 
                       stringsAsFactors = FALSE)

# Transform missing values
snp_data[(snp_data == "..") | (snp_data == "./.") | (snp_data == ".")] <- NA

# Reorder columns to match sample table
snp_data <- snp_data[, sample_ids]

# Load methylation data --------------------------------------------------------
cat("Loading methylation data...\n")
meth_data <- read.table(methylation_file, 
                        sep = '\t', 
                        header = TRUE, 
                        stringsAsFactors = FALSE)

# Set loci ID as rownames (assuming chr and pos columns exist)
if ("chr" %in% colnames(meth_data) & "pos" %in% colnames(meth_data)) {
  rownames(meth_data) <- paste0(meth_data$chr, "_", meth_data$pos)
}

################################################################################
# PART 2: PREPARE GENETIC DATA (SNPs)
################################################################################

cat("Processing SNP data...\n")

# Transpose SNP data (samples as rows, loci as columns)
snp_data_t <- data.frame(t(snp_data))
colnames(snp_data_t) <- gsub("X", "", colnames(snp_data_t))

# Build genind object
snp_genind <- df2genind(snp_data_t, 
                        ploidy = 2, 
                        ind.names = rownames(sample_table), 
                        pop = sample_table$population,  # Adjust column name
                        sep = "/")

cat("Genind object created:\n")
print(snp_genind)

# Scale allele frequencies -----------------------------------------------------
# Centers to mean zero and scales; replaces NAs with mean allele frequencies
scaled_snp <- scaleGen(snp_genind, NA.method = "mean")

cat("Scaled SNP matrix dimensions:", dim(scaled_snp), "\n")

################################################################################
# PART 3: PRINCIPAL COMPONENT ANALYSIS (PCA) ON GENETIC DATA
################################################################################

cat("Running PCA on genetic data...\n")

# Perform PCA
pca_snp <- dudi.pca(scaled_snp, 
                    cent = FALSE, 
                    scale = FALSE, 
                    scannf = FALSE)

# Plot variance explained by each PC
pdf(file.path(output_dir, "pca_variance_explained.pdf"), width = 10, height = 6)
barplot(100 * pca_snp$eig / sum(pca_snp$eig), 
        main = "PCA Eigenvalues - Variance Explained",
        col = heat.colors(50), 
        ylim = c(0, max(100 * pca_snp$eig / sum(pca_snp$eig)) + 1),
        ylab = "Percent of Genetic Variance Explained",
        xlab = "Principal Components")
dev.off()

# Calculate variance explained by first PCs
var_explained <- 100 * pca_snp$eig / sum(pca_snp$eig)
cat("Variance explained by PC1:", round(var_explained[1], 2), "%\n")
cat("Variance explained by PC2:", round(var_explained[2], 2), "%\n")
cat("Variance explained by PC3:", round(var_explained[3], 2), "%\n")
cat("Cumulative variance (PC1-3):", round(sum(var_explained[1:3]), 2), "%\n")

# Determine number of PCs to retain (e.g., explaining ~30% variance)
n_pcs <- 12  # Adjust based on your data
cat("Retaining", n_pcs, "PCs for RDA analysis\n")

# Re-run PCA with specified number of PCs
pca_snp <- dudi.pca(scaled_snp, 
                    cent = FALSE, 
                    scale = FALSE, 
                    scannf = FALSE, 
                    nf = n_pcs)

pcs <- round(pca_snp$li, 4)

# Plot PCA biplot
pdf(file.path(output_dir, "pca_biplot.pdf"), width = 10, height = 8)
col_palette <- rainbow(length(unique(sample_table$population)))
s.class(pca_snp$li, 
        pop(snp_genind),
        xax = 1, 
        yax = 2, 
        col = transp(col_palette, 0.6), 
        axesell = TRUE,
        addaxes = TRUE, 
        cstar = 0, 
        cpoint = 3, 
        grid = FALSE)
title(ylab = paste0("PC1 (", round(var_explained[1], 2), "%)"), line = 2)
title(xlab = paste0("PC2 (", round(var_explained[2], 2), "%)"), line = 1)
dev.off()

################################################################################
# PART 4: PREPARE METHYLATION DATA FOR RDA
################################################################################

cat("Processing methylation data...\n")

# Remove rows with NAs
info_cols <- 1:3  # Adjust based on your data structure (chr, pos, strand, etc.)
data_cols <- setdiff(1:ncol(meth_data), info_cols)

rows_with_na <- apply(meth_data[, data_cols], 1, function(x) any(is.na(x)))
meth_data_filt <- meth_data[!rows_with_na, ]

cat("Methylation data: removed", sum(rows_with_na), "rows with missing values\n")
cat("Remaining loci:", nrow(meth_data_filt), "\n")

# Transpose methylation data (samples as rows, loci as columns)
rda_meth <- data.frame(t(meth_data_filt[, data_cols]))
colnames(rda_meth) <- gsub("X", "", colnames(rda_meth))

################################################################################
# PART 5: REDUNDANCY ANALYSIS (RDA)
################################################################################

cat("Running RDA analysis...\n")

# Model 1: Population as predictor, conditioned on genetic PCs ----------------
rda_model1 <- rda(as.matrix(rda_meth) ~ population + Condition(as.matrix(pcs)),
                  data = sample_table)

cat("\n=== RDA Model 1: Population (conditioned on genetics) ===\n")
print(rda_model1)

# Perform permutation test
set.seed(123)  # For reproducibility
rda_anova1 <- anova(rda_model1, permutations = 999)
cat("\nPermutation test results:\n")
print(rda_anova1)

# Calculate adjusted R-squared
rsq1 <- RsquareAdj(rda_model1)
cat("\nR-squared:", round(rsq1$r.squared, 4), "\n")
cat("Adjusted R-squared:", round(rsq1$adj.r.squared, 4), "\n")

# Model 2: Treatment as predictor (if applicable) -----------------------------
if ("treatment" %in% colnames(sample_table)) {
  cat("\n=== RDA Model 2: Treatment (conditioned on genetics) ===\n")
  
  rda_model2 <- rda(as.matrix(rda_meth) ~ treatment + Condition(as.matrix(pcs)),
                    data = sample_table)
  
  print(rda_model2)
  
  set.seed(123)
  rda_anova2 <- anova(rda_model2, permutations = 999)
  print(rda_anova2)
  
  rsq2 <- RsquareAdj(rda_model2)
  cat("\nR-squared:", round(rsq2$r.squared, 4), "\n")
  cat("Adjusted R-squared:", round(rsq2$adj.r.squared, 4), "\n")
}

# Model 3: Interaction model (if applicable) ----------------------------------
if (all(c("country", "treatment") %in% colnames(sample_table))) {
  cat("\n=== RDA Model 3: Country + Treatment + Interaction ===\n")
  
  rda_model3 <- rda(as.matrix(rda_meth) ~ country + treatment + country:treatment + 
                      Condition(as.matrix(pcs)),
                    data = sample_table)
  
  print(rda_model3)
  
  set.seed(123)
  rda_anova3 <- anova(rda_model3, permutations = 999)
  print(rda_anova3)
  
  rsq3 <- RsquareAdj(rda_model3)
  cat("\nR-squared:", round(rsq3$r.squared, 4), "\n")
  cat("Adjusted R-squared:", round(rsq3$adj.r.squared, 4), "\n")
}

################################################################################
# PART 6: VARIANCE PARTITIONING
################################################################################

cat("\nCalculating variance explained by RDA axes...\n")

# Get variance explained by each constrained axis
rda_summary <- summary(eigenvals(rda_model1))
total_var <- sum(eigenvals(rda_model1))
var_per_axis <- eigenvals(rda_model1) / total_var

cat("Variance explained by first 5 RDA axes:\n")
print(round(var_per_axis[1:min(5, length(var_per_axis))], 4))

################################################################################
# PART 7: VISUALIZATION
################################################################################

cat("Creating RDA plots...\n")

# Extract site scores (sample coordinates in RDA space)
rda_sites <- data.frame(scores(rda_model1, choices = 1:4)$sites)
rda_sites <- cbind(rda_sites, sample_table)

# Define color scheme (adjust based on your grouping variables)
if ("treatment" %in% colnames(sample_table)) {
  color_var <- "treatment"
  color_palette <- c("C" = "#9DBF9E", "Cd" = "#FCB97D", "Cu" = "#FF5E5B")
} else {
  color_var <- "population"
  color_palette <- rainbow(length(unique(sample_table$population)))
}

# RDA Plot: Axis 1 vs Axis 2
p1 <- ggplot(rda_sites, aes(x = RDA1, y = RDA2)) +
  geom_point(aes(color = .data[[color_var]], 
                 shape = population),
             size = 4, 
             alpha = 0.7) +
  scale_color_manual(values = color_palette) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  xlab(paste0("RDA1 (", round(var_per_axis[1] * 100, 2), "%)")) +
  ylab(paste0("RDA2 (", round(var_per_axis[2] * 100, 2), "%)")) +
  coord_fixed() +
  ggtitle("RDA: Epigenetic Variation")

ggsave(file.path(output_dir, "rda_plot_axis1_vs_axis2.pdf"), 
       plot = p1, 
       width = 10, 
       height = 8)

# RDA Plot: Axis 1 vs Axis 3
if (ncol(rda_sites) >= 5) {
  p2 <- ggplot(rda_sites, aes(x = RDA1, y = RDA3)) +
    geom_point(aes(color = .data[[color_var]], 
                   shape = population),
               size = 4, 
               alpha = 0.7) +
    scale_color_manual(values = color_palette) +
    theme_bw(base_size = 14) +
    theme(legend.position = "right",
          legend.title = element_text(face = "bold")) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    xlab(paste0("RDA1 (", round(var_per_axis[1] * 100, 2), "%)")) +
    ylab(paste0("RDA3 (", round(var_per_axis[3] * 100, 2), "%)")) +
    coord_fixed() +
    ggtitle("RDA: Epigenetic Variation")
  
  ggsave(file.path(output_dir, "rda_plot_axis1_vs_axis3.pdf"), 
         plot = p2, 
         width = 10, 
         height = 8)
}

################################################################################
# PART 8: SAVE RESULTS
################################################################################

cat("Saving results...\n")

# Save RDA model objects
saveRDS(rda_model1, file.path(output_dir, "rda_model_population.rds"))

# Save site scores
write.table(rda_sites, 
            file.path(output_dir, "rda_site_scores.txt"),
            sep = "\t", 
            quote = FALSE, 
            row.names = TRUE)

# Save variance explained
variance_table <- data.frame(
  Axis = paste0("RDA", 1:length(var_per_axis)),
  Variance = var_per_axis,
  Percentage = var_per_axis * 100
)

write.table(variance_table, 
            file.path(output_dir, "rda_variance_explained.txt"),
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

cat("\nAnalysis complete! Results saved in:", output_dir, "\n")

################################################################################
# END OF SCRIPT
################################################################################