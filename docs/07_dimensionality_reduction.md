# (PART) DIMENSIONALITY REDUCTION {-}

# Dimensionality Reduction Techniques
Dimensionality reduction is a crucial step in microbiome analysis, aiming to simplify datasets while retaining essential information. This section explores techniques for reducing the number of features or variables in a phyloseq object, enabling more efficient analysis and visualization of microbial community data.

We delve into commonly used dimensionality reduction methods such as PCA, MDS, and t-SNE. These techniques utilize distance matrices obtained from the previous step to visualize and explore the structure of microbiome data in lower-dimensional spaces. The subsection covers the technical details of each method and provides examples of their application in microbiome data analysis.




## Prerequisite {-}
- Preprocessed R objects
- Preprocessed phyloseq objects

> Either of the R or phyloseq objects can be used as input data.

## Import libraries and data objects {-}

```r
set.seed(110912)

library(tidyverse, suppressPackageStartupMessages())
library(phyloseq)
library(microbiome)
library(vegan)

cat("\nSaved RData objects\n\n")

Saved RData objects
load("../imap-data-preparation/data/phyloseq_raw_rel_psextra_df_objects.rda", verbose = T)
Loading objects:
  ps_raw
  ps_rel
  psextra_raw
  psextra_rel
  ps_df
load("../imap-data-preparation/data/ps_transformed.rda", verbose = T)
Loading objects:
  ps_asin
  ps_identity
  ps_compositional
  ps_z_otu
  ps_z_sample
  ps_log10
  ps_log10p
  ps_clr
  ps_shift
  ps_scale
load("../imap-data-preparation/data/bray_distances.rda", verbose = T)
Loading objects:
  ps_asin_bray_dist
  ps_compositional_bray_dist
  ps_z_otu_bray_dist
  ps_z_sample_bray_dist
  ps_log10_bray_dist
  ps_log10p_bray_dist
  ps_clr_bray_dist
  ps_shift_bray_dist
  ps_scale_bray_dist
```

## PCA (Principal Component Analysis)

PCA is a dimensionality reduction technique commonly used to explore patterns and structure in high-dimensional data. It works by transforming the original variables into a new set of uncorrelated variables called principal components. These principal components capture the maximum variance in the data, allowing for visualization and exploration of the major patterns or clusters within the dataset.


```r
# Sample R code for PCA with distance matrices
library(stats)

# Assuming 'distance_matrix' is your distance matrix between samples
pca_asin_bray_metrics <- prcomp(ps_asin_bray_dist)
```


## MDS (Multidimensional Scaling)

MDS is a method for visualizing the similarities or dissimilarities between samples in a lower-dimensional space. It aims to preserve the pairwise distances or dissimilarities between samples as much as possible, representing the data in a reduced-dimensional space while maintaining the overall structure. MDS is particularly useful for exploring the overall structure of complex datasets and identifying sample clusters or gradients based on ecological distances or dissimilarity measures.


```r
# Sample R code for MDS with distance matrices
library(stats)

# Assuming 'distance_matrix' is your distance matrix between samples
mds_asin_bray_metrics <- cmdscale(ps_asin_bray_dist)
```

## PCoA (Principal Coordinates Analysis)
Principal Coordinates Analysis (PCoA), also known as classical multidimensional scaling (MDS), is a dimensionality reduction technique commonly used in ecology and other fields to visualize and explore patterns in distance or dissimilarity data. PCoA takes a distance or dissimilarity matrix as input and transforms it into a set of orthogonal axes (principal coordinates) that capture the maximum variation in the data. These principal coordinates represent a lower-dimensional representation of the original data while preserving the pairwise distances between samples as much as possible. PCoA is particularly useful for visualizing and interpreting relationships between samples based on their similarities or dissimilarities, facilitating the identification of clusters, gradients, or other patterns in complex datasets


```r
# Sample R code for MDS with distance matrices
library(vegan)

# Assuming 'distance_matrix' is your distance matrix between samples
pcoa_asin_bray_metrics <- capscale(ps_asin_bray_dist ~ 1)
```



## t-SNE (t-Distributed Stochastic Neighbor Embedding)

t-SNE is a nonlinear dimensionality reduction technique that aims to preserve both local and global relationships between samples in a lower-dimensional space. It is particularly effective for visualizing high-dimensional data in two or three dimensions, capturing fine-scale structure or clustering within the dataset. t-SNE works by modeling the similarities between samples in the high-dimensional space and mapping them to a lower-dimensional space while minimizing the Kullback-Leibler divergence between the two distributions. It is commonly used for exploratory data analysis and visualization of complex datasets, including microbiome abundance data.


```r
# Sample R code for t-SNE with distance matrices
library(Rtsne)

# Assuming 'distance_matrix' is your distance matrix between samples
tsne_asin_bray_metrics <- Rtsne(as.dist(ps_asin_bray_dist))
```



```r
save(pca_asin_bray_metrics, mds_asin_bray_metrics, pcoa_asin_bray_metrics, tsne_asin_bray_metrics, file = "data/reduced_dimension.rda")
```


## Non-negative Matrix Factorization (NMF)

NMF is a dimensionality reduction technique that factorizes a non-negative matrix into two lower-dimensional matrices with non-negative elements. In microbiome analysis, NMF can be applied to identify latent features or microbial signatures within complex microbial community data, facilitating the identification of microbial co-occurrence patterns or functional groups.




## Redundancy Analysis (RDA)

RDA is a multivariate statistical technique used to explore the relationships between multiple response variables (e.g., microbial taxa) and explanatory variables (e.g., environmental or metadata factors). In microbiome analysis, RDA can be applied to identify microbial taxa associated with specific environmental gradients or experimental conditions, providing insights into the ecological drivers shaping microbial composition.


```r
# Sample R code for RDA with distance matrices
library(vegan)

# Assuming 'distance_matrix' is your distance matrix between samples
rda_metrics <- rda(ps_asin)
```


## Canonical Correspondence Analysis (CCA)

CCA is a multivariate technique similar to RDA but specifically designed for analyzing relationships between species data (e.g., microbial taxa) and environmental variables. In microbiome analysis, CCA can be used to identify microbial taxa that are most strongly correlated with environmental or metadata factors, helping to understand how microbial communities respond to different ecological conditions.




## Random Forests

Random Forests is a machine learning algorithm that can be applied for feature selection and dimensionality reduction in microbiome data. It works by constructing an ensemble of decision trees and aggregating their predictions to identify the most important features (e.g., microbial taxa) for predicting a response variable of interest (e.g., disease status). Random Forests can be useful for identifying key microbial biomarkers associated with specific outcomes or conditions.




## Deep Learning

Deep learning techniques, such as deep neural networks, can also be applied for dimensionality reduction and feature extraction in microbiome data. These techniques have the ability to learn complex patterns and representations directly from the data, potentially uncovering novel insights into microbial community structure and function. However, deep learning methods may require large amounts of data and computational resources, and their interpretability can be challenging.



These methods provide alternative approaches for reducing the dimensionality of microbiome data and extracting meaningful information from complex microbial community datasets. Depending on the specific research question and characteristics of the data, different methods may be more suitable or informative.


# Taxa Prevalence and Detection

- Taxa prevalence and detection analysis quantify how often a taxon is present or detected across the samples under study.
- High prevalence indicates that the taxon is commonly found across many samples.
- Low prevalence suggests that the taxon is more sporadic or rare.
- Analyzing taxa prevalence can provide insights into the distribution of different microbial taxa and their potential ecological roles.

## Understand sample and taxa sizes

```r
library(phyloseq)

# Min and max sample and taxa sums: These values are useful when setting up some threshholds.

cat("Minimum sample sums:", min(data.frame(sample_sums(ps_raw))), "\n\n")
Minimum sample sums: 1776 
cat("Maximum sample sums:", max(data.frame(sample_sums(ps_raw))), "\n\n")
Maximum sample sums: 28883 
cat("Minimum taxa sums:", min(data.frame(taxa_sums(ps_raw))), "\n\n")
Minimum taxa sums: 0 
cat("Maximum taxa sums:", max(data.frame(taxa_sums(ps_raw))), "\n\n")
Maximum taxa sums: 873314 
```


## Subsetting Microbiome Data by Taxa Prevalence


```r
library(phyloseq)
library(microbiome)

# Define a function to subset microbiome data based on taxa prevalence
subset_by_prevalence <- function(ps_raw, prevalence_threshold) {
  # Calculate taxa prevalence across samples
  taxa_prevalence <- prevalence(ps_raw)
  
  # Subset taxa based on prevalence threshold
  taxa_to_keep <- rownames(taxa_prevalence)[taxa_prevalence >= prevalence_threshold]
  
  # Subset microbiome data
  ps_subset <- prune_taxa(taxa_to_keep, ps_raw)
  
  return(ps_subset)
}

# Example usage: Subset microbiome data with a prevalence threshold of 0.1 (10%)
ps_subset <- subset_by_prevalence(ps_raw, 0.9)
```

