# (PART) DATA TRANSFORMATION {-}

# Data Tranformation from Phyloseq Objects {#data-tarnsformation}

The process of data transformation in microbiome analysis involves converting raw or relative abundance values into matrices that are suitable for further analysis. Data transformation focuses on modifying the data itself. In this context, we will explore various methods used for this transformation on a phyloseq object.





## Import processed phyloseq object

```r
load("../imap-data-preparation/data/external_ps_objects.rda", verbose = TRUE) 
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
                                    Sample-1 Sample-2 Sample-3 Sample-4
Clostridium ramosum et rel.                2        2        2        2
Staphylococcus                             0        0        0        0
Papillibacter cinnamivorans et rel.        9       59       55       45
Campylobacter                              4        4        4        4
Micrococcaceae                             0        0        0        0
Prevotella ruminicola et rel.              1        2        3        2
Bryantella formatexigens et rel.          78      106      126       47
Veillonella                                8       19        6        5
Catenibacterium mitsuokai et rel.          0        5        3        2
Uncultured Clostridiales II               27      205      284      155
                                    Sample-5 Sample-6 Sample-7 Sample-8
Clostridium ramosum et rel.                2        2        3        2
Staphylococcus                             0        0        0        0
Papillibacter cinnamivorans et rel.      108      110       58       55
Campylobacter                              4        4        4        4
Micrococcaceae                             0        0        0        0
Prevotella ruminicola et rel.              9        4        2        8
Bryantella formatexigens et rel.          54       52       73       43
Veillonella                                7        8       11       11
Catenibacterium mitsuokai et rel.          0        8        7        1
Uncultured Clostridiales II               30       35      175       45
                                    Sample-9 Sample-10
Clostridium ramosum et rel.                2         2
Staphylococcus                             0         0
Papillibacter cinnamivorans et rel.       99        20
Campylobacter                              4         4
Micrococcaceae                             0         0
Prevotella ruminicola et rel.              1         2
Bryantella formatexigens et rel.          41        48
Veillonella                               71        74
Catenibacterium mitsuokai et rel.          2         1
Uncultured Clostridiales II               51       141
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
Clostridium ramosum et rel.                2        2        2
Staphylococcus                             0        0        0
Papillibacter cinnamivorans et rel.        9       59       55
Campylobacter                              4        4        4
Micrococcaceae                             0        0        0
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
Clostridium ramosum et rel.         0.0002364346 9.484966e-05 6.924488e-05
Staphylococcus                      0.0000000000 0.000000e+00 0.000000e+00
Papillibacter cinnamivorans et rel. 0.0010639556 2.798065e-03 1.904234e-03
Campylobacter                       0.0004728691 1.896993e-04 1.384898e-04
Micrococcaceae                      0.0000000000 0.000000e+00 0.000000e+00
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
Clostridium ramosum et rel.         0.017525 0.011099 0.009484 0.014005
Staphylococcus                      0.000000 0.000000 0.000000 0.000000
Papillibacter cinnamivorans et rel. 0.037182 0.060320 0.049752 0.066477
Campylobacter                       0.024785 0.015697 0.013412 0.019806
Micrococcaceae                      0.000000 0.000000 0.000000 0.000000
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
Clostridium ramosum et rel.         0.0002364346 9.484966e-05 6.924488e-05
Staphylococcus                      0.0000000000 0.000000e+00 0.000000e+00
Papillibacter cinnamivorans et rel. 0.0010639556 2.798065e-03 1.904234e-03
Campylobacter                       0.0004728691 1.896993e-04 1.384898e-04
Micrococcaceae                      0.0000000000 0.000000e+00 0.000000e+00
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
                                        Sample-1     Sample-2     Sample-3
Clostridium ramosum et rel.         -0.584998298 -0.584998298 -0.584998298
Staphylococcus                      -0.268583981 -0.268583981 -0.268583981
Papillibacter cinnamivorans et rel. -1.141005100  1.316483388  1.221856185
Campylobacter                       -0.008355085 -0.008355085 -0.008355085
Micrococcaceae                       0.000000000  0.000000000  0.000000000
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
                                       Sample-1   Sample-2   Sample-3
Clostridium ramosum et rel.         -0.56599933 -0.7027117 -0.6855514
Staphylococcus                      -1.15179151 -1.2167251 -1.1804675
Papillibacter cinnamivorans et rel.  0.07597223  0.6989171  0.6329215
Campylobacter                       -0.29362148 -0.4637090 -0.4554285
Micrococcaceae                      -1.15179151 -1.2167251 -1.1804675
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
                                     Sample-1  Sample-2  Sample-3
Clostridium ramosum et rel.         0.4771213 0.4771213 0.4771213
Staphylococcus                      0.0000000 0.0000000 0.0000000
Papillibacter cinnamivorans et rel. 1.0000000 1.7781513 1.7481880
Campylobacter                       0.6989700 0.6989700 0.6989700
Micrococcaceae                      0.0000000 0.0000000 0.0000000
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
                                     Sample-1  Sample-2  Sample-3
Clostridium ramosum et rel.         0.4771213 0.4771213 0.4771213
Staphylococcus                      0.0000000 0.0000000 0.0000000
Papillibacter cinnamivorans et rel. 1.0000000 1.7781513 1.7481880
Campylobacter                       0.6989700 0.6989700 0.6989700
Micrococcaceae                      0.0000000 0.0000000 0.0000000
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
                                      Sample-1   Sample-2   Sample-3
Clostridium ramosum et rel.         -0.8704964 -1.4908696 -1.5055455
Staphylococcus                      -3.5554744 -3.3594489 -3.1149834
Papillibacter cinnamivorans et rel.  0.5790592  1.7320473  1.5945468
Campylobacter                       -0.2120560 -0.8780341 -0.9177588
Micrococcaceae                      -3.5554744 -3.3594489 -3.1149834
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
Clostridium ramosum et rel.                3        3        3
Staphylococcus                             1        1        1
Papillibacter cinnamivorans et rel.       10       60       56
Campylobacter                              5        5        5
Micrococcaceae                             1        1        1
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
Clostridium ramosum et rel.                2        2        2
Staphylococcus                             0        0        0
Papillibacter cinnamivorans et rel.        9       59       55
Campylobacter                              4        4        4
Micrococcaceae                             0        0        0
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


