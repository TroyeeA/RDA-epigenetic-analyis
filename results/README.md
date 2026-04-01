# RDA Epigenetic Analysis

Redundancy Analysis (RDA) pipeline for analyzing genotype-environment associations using SNP and DNA methylation data.

## Overview

This repository contains R scripts for performing Redundancy Analysis (RDA) to investigate how genetic and epigenetic variation respond to environmental predictors. The analysis includes:

- Principal Component Analysis (PCA) on SNP data
- RDA on methylation data conditioned on genetic structure
- Variance partitioning
- Visualization of results

## Requirements

### R Packages
```r
install.packages(c(
  "data.table",
  "ggplot2",
  "tidyverse",
  "gridExtra",
  "hierfstat",
  "impute",
  "psych",
  "vegan"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("adegenet")
```

## Input Data Format

### 1. Sample Metadata (`sample_design.txt`)
#Tab-delimited file with samples as rows:
sample_id population treatment country 

### 2. SNP Data (`snp_data.txt`)
#Tab-delimited file with loci as rows and samples as columns. Genotypes in format `A/T`, `./. ` for missing.

### 3. Methylation Data (`methylation_data.bed`)
BED-format file with columns: chr, pos, strand, followed by methylation values for each sample.

## Usage

### Basic Analysis
```r
# 1. Clone this repository
# 2. Place your data files in the data/ directory
# 3. Modify input paths in the script
# 4. Run the analysis

source("scripts/rda_analysis.R")
```

### Customize Analysis

Edit these parameters in the script:
```r
# Number of PCs to retain
n_pcs <- 12

# Predictor variable(s)
# Change 'population' to your predictor of interest
rda_model <- rda(methylation ~ population + Condition(PCs), data = metadata)
```

## Output

Results are saved in the `results/` directory:

- `pca_variance_explained.pdf` - Scree plot of PCA
- `pca_biplot.pdf` - PCA visualization
- `rda_plot_axis1_vs_axis2.pdf` - RDA ordination plot
- `rda_site_scores.txt` - Sample coordinates in RDA space
- `rda_variance_explained.txt` - Variance explained by each RDA axis
- `rda_model_population.rds` - Saved RDA model object

## Theory

RDA is a multivariate linear regression where each response variable is regressed against predictors, followed by PCA of fitted values to produce canonical axes. These axes represent linear combinations of predictors that best explain variation in the response data.

**Key assumption**: Linear relationship between response (genotypes/methylation) and explanatory variables (environment).

## References

- Legendre P, Legendre L (2012) Numerical Ecology, 3rd edition. Elsevier, Amsterdam.
- Borcard D, Gillet F, Legendre P (2011) Numerical Ecology with R. Springer, New York.
- [RDA Tutorial](https://popgen.nescent.org/2018-03-27_RDA_GEA.html)

## Citation

If you use this code, please cite:
Anupoma N. Troyee (2024). RDA Epigenetic Analysis Pipeline. GitHub repository: [https://github.com/TroyeeA/RDA-epigenetic-analyis]


## License

MIT License - see LICENSE file

## Contact

For questions or issues, please open an issue on GitHub or contact: niloyatroyee@example.com

## Acknowledgments

- Original workflow adapted from [PopGen RDA tutorial](https://popgen.nescent.org/2018-03-27_RDA_GEA.html)
