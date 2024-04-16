# (PART) DISTANCES MATRICES {-}

# Distance Calculation in Microbiome Data Analysis
Calculating distance matrices between samples using appropriate metrics is a fundamental step in data preprocessing. This practice efficiently captures the similarity or dissimilarity between samples, laying the groundwork for subsequent analysis and dimensionality reduction techniques such as PCA, MDS, or t-SNE.




## Import processed and transformed data objects

```r
set.seed(110912)

library(phyloseq)
library(microbiome)
library(vegan)

load("data/ps_transformed.rda", verbose = TRUE)
```

```
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
```

```r
load("data/phyloseq_objects.rda", verbose = TRUE)
```

```
Loading objects:
  ps_GlobalPatterns
  ps_dietswap
  ps_caporaso
  ps_kostic_crc
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


## Confirm computed Bray-Curtis distance objects

```r
load("data/bray_distances.rda", verbose = T)
```

```
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

## Confirm computed Jaccard distance objects

```r
load("data/jaccard_distances.rda", verbose = T)
```

```
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


## UniFrac Distance (Weighted and Unweighted)
UniFrac distances quantify the dissimilarity between microbial communities based on the evolutionary distances along a phylogenetic tree. Unweighted UniFrac considers only the presence or absence of taxa, while weighted UniFrac also takes into account the relative abundance of taxa.



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


## Aitchison Distance
The Aitchison distance, also known as the compositional distance or log-ratio distance, is specifically designed for compositional data such as microbiome abundance data. It calculates the dissimilarity between two samples based on the log-ratio transformations of their compositional vectors, providing meaningful comparisons between samples while addressing the compositional constraints inherent in relative abundance data.

- Computing aitchison using `microViz` packages


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


