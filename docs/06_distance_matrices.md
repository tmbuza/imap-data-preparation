# (PART) DISTANCES MATRICES {-}

# Distance Calculation in Microbiome Data Analysis
Calculating distance matrices between samples using appropriate metrics is a fundamental step in data preprocessing. This practice efficiently captures the similarity or dissimilarity between samples, laying the groundwork for subsequent analysis and dimensionality reduction techniques such as PCA, MDS, or t-SNE.




## Import processed and transformed data objects

```r
set.seed(110912)

library(phyloseq)
library(microbiome)
library(vegan)

load("data/phyloseq_objects.rda", verbose = TRUE)
Loading objects:
  ps_GlobalPatterns
  ps_dietswap
  ps_caporaso
  ps_kostic_crc
load("data/ps_transformed.rda", verbose = TRUE)
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
load("data/bray_distances.rda", verbose = TRUE)
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
load("data/phyloseq_raw_rel_psextra_df_objects.rda", verbose = TRUE)
Loading objects:
  ps_raw
  ps_rel
  psextra_raw
  psextra_rel
  ps_df
```



## Bray-Curtis Dissimilarity
The Bray-Curtis dissimilarity calculates the dissimilarity between two samples based on the relative abundances of taxa. It considers both presence and abundance information but is not sensitive to differences in abundance magnitude.


```r
library(vegan)
library(dplyr)

ps_asin_bray_dist  <- vegdist(t(otu_table(ps_asin)), method="bray") %>%
  as.matrix()
ps_identity_bray_dist  <- vegdist(t(otu_table(ps_identity)), method="bray") %>%
  as.matrix 
ps_compositional_bray_dist  <- vegdist(t(otu_table(ps_compositional)), method="bray") %>%
  as.matrix 
ps_z_otu_bray_dist  <- vegdist(t(otu_table(ps_z_otu)), method="bray") %>%
  as.matrix 
ps_z_sample_bray_dist  <- vegdist(t(otu_table(ps_z_sample)), method="bray") %>%
  as.matrix 
ps_log10_bray_dist  <- vegdist(t(otu_table(ps_log10)), method="bray") %>%
  as.matrix 
ps_log10p_bray_dist  <- vegdist(t(otu_table(ps_log10p)), method="bray") %>%
  as.matrix 
ps_clr_bray_dist  <- vegdist(t(otu_table(ps_clr)), method="bray") %>%
  as.matrix 
ps_shift_bray_dist  <- vegdist(t(otu_table(ps_shift)), method="bray") %>%
  as.matrix 
ps_scale_bray_dist  <- vegdist(t(otu_table(ps_scale)), method="bray") %>%
  as.matrix 

save(ps_asin_bray_dist, 
     ps_compositional_bray_dist, 
     ps_z_otu_bray_dist, 
     ps_z_sample_bray_dist, 
     ps_log10_bray_dist, 
     ps_log10p_bray_dist, 
     ps_clr_bray_dist, 
     ps_shift_bray_dist, 
     ps_scale_bray_dist, 
     file = "data/bray_distances.rda")
```


## Jaccard Distance

The Jaccard distance measures the dissimilarity between two samples based on the presence or absence of taxa, ignoring their abundance information. It is calculated as the ratio of the number of taxa found in only one of the samples to the total number of taxa found in both samples.


```r
library(vegan)
library(dplyr)

ps_asin_jaccard_dist  <- vegdist(t(otu_table(ps_asin)), method="jaccard") %>%
  as.matrix()
ps_identity_jaccard_dist  <- vegdist(t(otu_table(ps_identity)), method="jaccard") %>%
  as.matrix 
ps_compositional_jaccard_dist  <- vegdist(t(otu_table(ps_compositional)), method="jaccard") %>%
  as.matrix 
ps_z_otu_jaccard_dist  <- vegdist(t(otu_table(ps_z_otu)), method="jaccard") %>%
  as.matrix 
ps_z_sample_jaccard_dist  <- vegdist(t(otu_table(ps_z_sample)), method="jaccard") %>%
  as.matrix 
ps_log10_jaccard_dist  <- vegdist(t(otu_table(ps_log10)), method="jaccard") %>%
  as.matrix 
ps_log10p_jaccard_dist  <- vegdist(t(otu_table(ps_log10p)), method="jaccard") %>%
  as.matrix 
ps_clr_jaccard_dist  <- vegdist(t(otu_table(ps_clr)), method="jaccard") %>%
  as.matrix 
ps_shift_jaccard_dist  <- vegdist(t(otu_table(ps_shift)), method="jaccard") %>%
  as.matrix 
ps_scale_jaccard_dist  <- vegdist(t(otu_table(ps_scale)), method="jaccard") %>%
  as.matrix 

save(ps_asin_jaccard_dist, 
     ps_compositional_jaccard_dist, 
     ps_z_otu_jaccard_dist, 
     ps_z_sample_jaccard_dist, 
     ps_log10_jaccard_dist, 
     ps_log10p_jaccard_dist, 
     ps_clr_jaccard_dist, 
     ps_shift_jaccard_dist, 
     ps_scale_jaccard_dist, 
     file = "data/jaccard_distances.rda")
```

Confirm computed Jaccard distance objects

```r
load("data/jaccard_distances.rda", verbose = T)
Loading objects:
  ps_asin_jaccard_dist
  ps_compositional_jaccard_dist
  ps_z_otu_jaccard_dist
  ps_z_sample_jaccard_dist
  ps_log10_jaccard_dist
  ps_log10p_jaccard_dist
  ps_clr_jaccard_dist
  ps_shift_jaccard_dist
  ps_scale_jaccard_dist
```



## Euclidean Distance
The Euclidean distance calculates the straight-line distance between two samples in a multidimensional space based on their abundance profiles. It considers differences in both presence and abundance of taxa.



```r
library(vegan)
library(dplyr)

ps_asin_euclidean_dist  <- vegdist(t(otu_table(ps_asin)), method="euclidean") %>%
  as.matrix()
ps_identity_euclidean_dist  <- vegdist(t(otu_table(ps_identity)), method="euclidean") %>%
  as.matrix 
ps_compositional_euclidean_dist  <- vegdist(t(otu_table(ps_compositional)), method="euclidean") %>%
  as.matrix 
ps_z_otu_euclidean_dist  <- vegdist(t(otu_table(ps_z_otu)), method="euclidean") %>%
  as.matrix 
ps_z_sample_euclidean_dist  <- vegdist(t(otu_table(ps_z_sample)), method="euclidean") %>%
  as.matrix 
ps_log10_euclidean_dist  <- vegdist(t(otu_table(ps_log10)), method="euclidean") %>%
  as.matrix 
ps_log10p_euclidean_dist  <- vegdist(t(otu_table(ps_log10p)), method="euclidean") %>%
  as.matrix 
ps_clr_euclidean_dist  <- vegdist(t(otu_table(ps_clr)), method="euclidean") %>%
  as.matrix 
ps_shift_euclidean_dist  <- vegdist(t(otu_table(ps_shift)), method="euclidean") %>%
  as.matrix 
ps_scale_euclidean_dist  <- vegdist(t(otu_table(ps_scale)), method="euclidean") %>%
  as.matrix 

save(ps_asin_euclidean_dist, 
     ps_compositional_euclidean_dist, 
     ps_z_otu_euclidean_dist, 
     ps_z_sample_euclidean_dist, 
     ps_log10_euclidean_dist, 
     ps_log10p_euclidean_dist, 
     ps_clr_euclidean_dist, 
     ps_shift_euclidean_dist, 
     ps_scale_euclidean_dist, 
     file = "data/euclidean_distances.rda")
```


Confirm computed Euclidean distance objects

```r
load("data/euclidean_distances.rda", verbose = T)
Loading objects:
  ps_asin_euclidean_dist
  ps_compositional_euclidean_dist
  ps_z_otu_euclidean_dist
  ps_z_sample_euclidean_dist
  ps_log10_euclidean_dist
  ps_log10p_euclidean_dist
  ps_clr_euclidean_dist
  ps_shift_euclidean_dist
  ps_scale_euclidean_dist
```


## Aitchison Distance
The Aitchison distance, also known as the compositional distance or log-ratio distance, is specifically designed for compositional data such as microbiome abundance data. It calculates the dissimilarity between two samples based on the log-ratio transformations of their compositional vectors, providing meaningful comparisons between samples while addressing the compositional constraints inherent in relative abundance data.



```r
library(microViz)
library(dplyr)

# ps_identity_aitchison_dists <- ps_asin %>% dist_calc("aitchison")
ps_GlobalPatterns_aitchison_dists <- ps_GlobalPatterns %>% dist_calc("aitchison")
ps_dietswap_aitchison_dists <- ps_dietswap %>% dist_calc("aitchison")
ps_kostic_crc_aitchison_dists <- ps_kostic_crc %>% dist_calc("aitchison")

save(ps_GlobalPatterns_aitchison_dists, 
     ps_dietswap_aitchison_dists, 
     ps_kostic_crc_aitchison_dists, 
     file = "data/aitchison_distances.rda")
```



Confirm computed Aitchison distance objects

```r
load("data/aitchison_distances.rda", verbose = T)
Loading objects:
  ps_GlobalPatterns_aitchison_dists
  ps_dietswap_aitchison_dists
  ps_kostic_crc_aitchison_dists
```


# What's Next: Utilizing Transformed Data {#utilizing-transformed-data}

After computing distances using various methods, such as Bray-Curtis, Jaccard, Euclidean, and Aitchison, we can harness these transformed data in downstream analyses. These analyses may encompass:

- **Beta Diversity Analysis:** Comparing microbial community compositions between samples to elucidate differences or similarities.
- **Multivariate Statistical Analysis:** Assessing the significance of differences between groups of samples using methods like PERMANOVA or MANOVA.
- **Data Visualization:** Identifying patterns or clusters within the microbial community data through visualization techniques such as PCoA or NMDS plots.
- **Correlation Analysis:** Exploring relationships between microbial taxa and environmental variables to uncover potential ecological associations.
- **Dimensionality Reduction:** Employing techniques like Principal Component Analysis (PCA), t-Distributed Stochastic Neighbor Embedding (t-SNE), or Uniform Manifold Approximation and Projection (UMAP) to visualize high-dimensional microbiome data in lower-dimensional spaces.

By leveraging the transformed data from distance calculations, we can gain deeper insights into microbial community dynamics and their relationships with various environmental factors or experimental conditions.


