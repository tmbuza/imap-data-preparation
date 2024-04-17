# (PART) DATA TRANSFORMATION {-}

# Data Tranformation from Phyloseq Objects {#data-tarnsformation}

The process of data transformation in microbiome analysis involves converting raw or relative abundance values into matrices that are suitable for further analysis. Data transformation focuses on modifying the data itself. In this context, we will explore various methods used for this transformation on a phyloseq object.





## Import processed phyloseq object

```r
load("../imap-data-preparation/data/phyloseq_objects.rda", verbose = TRUE) 
```

```
Loading objects:
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
Staphylococcus                          0        0        0        0        0
Uncultured Clostridiales I             19      120     1183      399       19
Alcaligenes faecalis et rel.            1        2        3        2        2
Bacteroides ovatus et rel.            427       15       31       46       29
Uncultured Bacteroidetes                2        2        1        5        1
Dialister                               5       23        6       24        6
Clostridium thermocellum et rel.        0        0        0        0        0
Novosphingobium                         0        0        0        0        0
Ruminococcus gnavus et rel.            15       56       65       41       47
Sutterella wadsworthia et rel.         77       11       44        7       55
                                 Sample-6 Sample-7 Sample-8 Sample-9 Sample-10
Staphylococcus                          0        0        0        0         0
Uncultured Clostridiales I             43      137       70       59       315
Alcaligenes faecalis et rel.            2        2        2        2         7
Bacteroides ovatus et rel.             13       15       19        7        11
Uncultured Bacteroidetes                0        1        1        0         2
Dialister                              38       37       25       45         4
Clostridium thermocellum et rel.        0        0        0        0         0
Novosphingobium                         0        0        0        0         0
Ruminococcus gnavus et rel.            61      101       19       21        19
Sutterella wadsworthia et rel.         14       13       10        6        10
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
Staphylococcus                      0        0        0
Uncultured Clostridiales I         19      120     1183
Alcaligenes faecalis et rel.        1        2        3
Bacteroides ovatus et rel.        427       15       31
Uncultured Bacteroidetes            2        2        1
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
Staphylococcus               0.0000000000 0.000000e+00 0.000000e+00
Uncultured Clostridiales I   0.0022461284 5.690980e-03 4.095835e-02
Alcaligenes faecalis et rel. 0.0001182173 9.484966e-05 1.038673e-04
Bacteroides ovatus et rel.   0.0504787800 7.113725e-04 1.073296e-03
Uncultured Bacteroidetes     0.0002364346 9.484966e-05 3.462244e-05
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
Staphylococcus               0.000000 0.000000 0.000000 0.000000
Uncultured Clostridiales I   0.054038 0.086080 0.232740 0.199116
Alcaligenes faecalis et rel. 0.012391 0.011099 0.011615 0.014005
Bacteroides ovatus et rel.   0.258934 0.030401 0.037345 0.067213
Uncultured Bacteroidetes     0.017525 0.011099 0.006706 0.022144
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
Staphylococcus               0.0000000000 0.000000e+00 0.000000e+00
Uncultured Clostridiales I   0.0022461284 5.690980e-03 4.095835e-02
Alcaligenes faecalis et rel. 0.0001182173 9.484966e-05 1.038673e-04
Bacteroides ovatus et rel.   0.0504787800 7.113725e-04 1.073296e-03
Uncultured Bacteroidetes     0.0002364346 9.484966e-05 3.462244e-05
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
                                Sample-1    Sample-2    Sample-3
Staphylococcus               -0.26858398 -0.26858398 -0.26858398
Uncultured Clostridiales I   -0.73136690  1.04729654  3.30105031
Alcaligenes faecalis et rel. -0.93076516 -0.44189484 -0.09503583
Bacteroides ovatus et rel.    1.80487522 -1.13754389 -0.51697242
Uncultured Bacteroidetes      0.07483888  0.07483888 -0.50384818
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
Staphylococcus               -1.1517915 -1.21672513 -1.1804675
Uncultured Clostridiales I    0.4455659  1.02710632  2.0075089
Alcaligenes faecalis et rel. -0.7821978 -0.70271169 -0.5559529
Bacteroides ovatus et rel.    2.0789997  0.08050037  0.3808190
Uncultured Bacteroidetes     -0.5659993 -0.70271169 -0.8682102
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
Staphylococcus               0.0000000 0.0000000 0.000000
Uncultured Clostridiales I   1.3010300 2.0827854 3.073352
Alcaligenes faecalis et rel. 0.3010300 0.4771213 0.602060
Bacteroides ovatus et rel.   2.6314438 1.2041200 1.505150
Uncultured Bacteroidetes     0.4771213 0.4771213 0.301030
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
Staphylococcus               0.0000000 0.0000000 0.000000
Uncultured Clostridiales I   1.3010300 2.0827854 3.073352
Alcaligenes faecalis et rel. 0.3010300 0.4771213 0.602060
Bacteroides ovatus et rel.   2.6314438 1.2041200 1.505150
Uncultured Bacteroidetes     0.4771213 0.4771213 0.301030
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
Staphylococcus               -3.5554744 -3.3594489 -3.114983
Uncultured Clostridiales I    1.3178115  2.4388711  4.654395
Alcaligenes faecalis et rel. -1.4976474 -1.4908696 -1.169073
Bacteroides ovatus et rel.    4.4228219  0.3804359  1.028151
Uncultured Bacteroidetes     -0.8704964 -1.4908696 -2.016371
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
Staphylococcus                      1        1        1
Uncultured Clostridiales I         20      121     1184
Alcaligenes faecalis et rel.        2        3        4
Bacteroides ovatus et rel.        428       16       32
Uncultured Bacteroidetes            3        3        2
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
Staphylococcus                      0        0        0
Uncultured Clostridiales I         19      120     1183
Alcaligenes faecalis et rel.        1        2        3
Bacteroides ovatus et rel.        427       15       31
Uncultured Bacteroidetes            2        2        1
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


