# (PART) DATA TRANSFORMATION {-}

# Data Tranformation from Phyloseq Objects {#data-tarnsformation}

The process of data transformation in microbiome analysis involves converting raw or relative abundance values into matrices that are suitable for further analysis. Data transformation focuses on modifying the data itself. In this context, we will explore various methods used for this transformation on a phyloseq object.





## Import processed phyloseq object

```r
load("../imap-data-preparation/data/external/external_ps_objects.rda", verbose = TRUE) 
```

```
Loading objects:
  df_GlobalPatterns
  df_dietswap
  df_caporaso
  df_kostic_crc
  ps_GlobalPatterns
  ps_dietswap
  ps_caporaso
  ps_kostic_crc
```

## Raw Abundance


```r
ps_raw <- ps_dietswap
otu_table(ps_raw)[1:10, 1:10]
```

```
OTU Table:          [10 taxa and 10 samples]
                     taxa are rows
                                  Sample-1 Sample-2 Sample-3 Sample-4 Sample-5
Clostridium felsineum et rel.            0        0        0        0        0
Streptococcus mitis et rel.             41       41       54       21       43
Catenibacterium mitsuokai et rel.        0        5        3        2        0
Pseudomonas                              0        0        0        0        0
Bacteroides ovatus et rel.             427       15       31       46       29
Aeromonas                                0        0        0        0        0
Fusobacteria                             4        5        5        5        5
Akkermansia                             18       97       67      256       21
Campylobacter                            4        4        4        4        4
Lactococcus                              1        4        2       81        1
                                  Sample-6 Sample-7 Sample-8 Sample-9 Sample-10
Clostridium felsineum et rel.            0        0        0        0         0
Streptococcus mitis et rel.            111       23       44      243        37
Catenibacterium mitsuokai et rel.        8        7        1        2         1
Pseudomonas                              0        0        0        0         0
Bacteroides ovatus et rel.              13       15       19        7        11
Aeromonas                                0        0        0        0         0
Fusobacteria                             8        5        5        4         5
Akkermansia                             16       26       30       19       125
Campylobacter                            4        4        4        4         4
Lactococcus                              6        1        1        7         1
```

##  No Transformation

> Similar to raw abundance


```r
(ps_identity <- microbiome::transform(ps_raw, 'identity'))
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
```

```r
cat("\n\n")
```

```r
otu_table(ps_identity)[1:5, 1:3]
```

```
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                                  Sample-1 Sample-2 Sample-3
Clostridium felsineum et rel.            0        0        0
Streptococcus mitis et rel.             41       41       54
Catenibacterium mitsuokai et rel.        0        5        3
Pseudomonas                              0        0        0
Bacteroides ovatus et rel.             427       15       31
```

## Relative abundance


```r
(ps_rel = phyloseq::transform_sample_counts(ps_raw, function(x){x / sum(x)}))
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
```

```r
otu_table(ps_rel)[1:5, 1:3]
```

```
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                                     Sample-1     Sample-2     Sample-3
Clostridium felsineum et rel.     0.000000000 0.0000000000 0.0000000000
Streptococcus mitis et rel.       0.004846909 0.0019444181 0.0018696119
Catenibacterium mitsuokai et rel. 0.000000000 0.0002371242 0.0001038673
Pseudomonas                       0.000000000 0.0000000000 0.0000000000
Bacteroides ovatus et rel.        0.050478780 0.0007113725 0.0010732957
```

## Arc sine (asin) transformation

- Typically used when dealing with proportional and percentages. 
- Proportionals range from 0 to 1
- Percentages range from 0 to 100
- The Metaphlan3 relative abundances are in percentages! That means the column totals in 100.


```r
x = otu_table(ps_rel)
y <- x/max(x)

ps_asin <- round(asin(sqrt(y)), 6)

ps_asin <- as.matrix(ps_asin)
ps_asin[1:5, 1:4]
```

```
OTU Table:          [5 taxa and 4 samples]
                     taxa are rows
                                  Sample-1 Sample-2 Sample-3 Sample-4
Clostridium felsineum et rel.     0.000000 0.000000 0.000000 0.000000
Streptococcus mitis et rel.       0.079426 0.050275 0.049297 0.045395
Catenibacterium mitsuokai et rel. 0.000000 0.017550 0.011615 0.014005
Pseudomonas                       0.000000 0.000000 0.000000 0.000000
Bacteroides ovatus et rel.        0.258934 0.030401 0.037345 0.067213
```


## Compositional Version
Compositional data represents relative proportions or percentages of different microbial taxa within a sample, rather than the absolute abundance of each taxon. This transformation is necessary because raw abundance data is typically affected by various factors, such as sequencing depth, which can lead to spurious correlations and biases in downstream analyses.


```r
(ps_compositional <- microbiome::transform(ps_raw, 'compositional'))
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
```

```r
cat("\n\n")
```

```r
otu_table(ps_compositional)[1:5, 1:3]
```

```
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                                     Sample-1     Sample-2     Sample-3
Clostridium felsineum et rel.     0.000000000 0.0000000000 0.0000000000
Streptococcus mitis et rel.       0.004846909 0.0019444181 0.0018696119
Catenibacterium mitsuokai et rel. 0.000000000 0.0002371242 0.0001038673
Pseudomonas                       0.000000000 0.0000000000 0.0000000000
Bacteroides ovatus et rel.        0.050478780 0.0007113725 0.0010732957
```

## Z-transform for OTUs


```r
(ps_z_otu <- microbiome::transform(ps_raw, 'Z', 'OTU'))
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
```

```r
cat("\n\n")
```

```r
otu_table(ps_z_otu)[1:5, 1:3]
```

```
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                                    Sample-1   Sample-2   Sample-3
Clostridium felsineum et rel.      0.0000000  0.0000000  0.0000000
Streptococcus mitis et rel.        0.5316419  0.5316419  0.7956010
Catenibacterium mitsuokai et rel. -1.0236796  0.9649155  0.5149076
Pseudomonas                       -0.2050931 -0.2050931 -0.2050931
Bacteroides ovatus et rel.         1.8048752 -1.1375439 -0.5169724
```

## Z-transform for Samples


```r
(ps_z_sample <- microbiome::transform(ps_raw, 'Z', 'sample'))
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
```

```r
cat("\n\n")
```

```r
otu_table(ps_z_sample)[1:5, 1:3]
```

```
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                                    Sample-1    Sample-2   Sample-3
Clostridium felsineum et rel.     -1.1517915 -1.21672513 -1.1804675
Streptococcus mitis et rel.        0.8411751  0.53203778  0.6248043
Catenibacterium mitsuokai et rel. -1.1517915 -0.37840531 -0.5559529
Pseudomonas                       -1.1517915 -1.21672513 -1.1804675
Bacteroides ovatus et rel.         2.0789997  0.08050037  0.3808190
```

## Log10 Transform


```r
(ps_log10 <- microbiome::transform(ps_raw, 'log10'))
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
```

```r
cat("\n\n")
```

```r
otu_table(ps_log10)[1:5, 1:3]
```

```
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                                  Sample-1  Sample-2 Sample-3
Clostridium felsineum et rel.     0.000000 0.0000000 0.000000
Streptococcus mitis et rel.       1.623249 1.6232493 1.740363
Catenibacterium mitsuokai et rel. 0.000000 0.7781513 0.602060
Pseudomonas                       0.000000 0.0000000 0.000000
Bacteroides ovatus et rel.        2.631444 1.2041200 1.505150
```

## Log10p Transform


```r
(ps_log10p <- microbiome::transform(ps_raw, 'log10p'))
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
```

```r
cat("\n\n")
```

```r
otu_table(ps_log10p)[1:5, 1:3]
```

```
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                                  Sample-1  Sample-2 Sample-3
Clostridium felsineum et rel.     0.000000 0.0000000 0.000000
Streptococcus mitis et rel.       1.623249 1.6232493 1.740363
Catenibacterium mitsuokai et rel. 0.000000 0.7781513 0.602060
Pseudomonas                       0.000000 0.0000000 0.000000
Bacteroides ovatus et rel.        2.631444 1.2041200 1.505150
```

## CLR Transform

- Note that small pseudocount is added if data contains zeroes

```r
(ps_clr <- microbiome::transform(ps_raw, 'clr'))
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
```

```r
cat("\n\n")
```

```r
otu_table(ps_clr)[1:5, 1:3]
```

```
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                                   Sample-1   Sample-2  Sample-3
Clostridium felsineum et rel.     -3.555474 -3.3594489 -3.114983
Streptococcus mitis et rel.        2.082832  1.3707778  1.576365
Catenibacterium mitsuokai et rel. -3.555474 -0.6717570 -1.169073
Pseudomonas                       -3.555474 -3.3594489 -3.114983
Bacteroides ovatus et rel.         4.422822  0.3804359  1.028151
```


## Shift the baseline


```r
(ps_shift <- microbiome::transform(ps_raw, 'shift', shift=1))
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
```

```r
cat("\n\n")
```

```r
otu_table(ps_shift)[1:5, 1:3]
```

```
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                                  Sample-1 Sample-2 Sample-3
Clostridium felsineum et rel.            1        1        1
Streptococcus mitis et rel.             42       42       55
Catenibacterium mitsuokai et rel.        1        6        4
Pseudomonas                              1        1        1
Bacteroides ovatus et rel.             428       16       32
```

## Data Scaling


```r
(ps_scale <- microbiome::transform(ps_raw, 'scale', scale=1))
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
```

```r
cat("\n\n")
```

```r
otu_table(ps_scale)[1:5, 1:3]
```

```
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                                  Sample-1 Sample-2 Sample-3
Clostridium felsineum et rel.            0        0        0
Streptococcus mitis et rel.             41       41       54
Catenibacterium mitsuokai et rel.        0        5        3
Pseudomonas                              0        0        0
Bacteroides ovatus et rel.             427       15       31
```

## Save transformed objects


```r
save(
  ps_asin, 
  ps_identity,
  ps_compositional, 
  ps_z_otu, 
  ps_z_sample, 
  ps_log10, 
  ps_log10p, 
  ps_clr, 
  ps_shift, 
  ps_scale, 
  file = "data/ps_transformed.rda")
```

## Confirm saved object


```r
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


