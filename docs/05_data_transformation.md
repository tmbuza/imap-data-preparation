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
Akkermansia                             18       97       67      256       21
Atopobium                                0        1        0        1        0
Megasphaera elsdenii et rel.             3       14       46       31        7
Oceanospirillum                          1        2        1        2        1
Prevotella tannerae et rel.             63       19       16       15       15
Catenibacterium mitsuokai et rel.        0        5        3        2        0
Helicobacter                             2        2        2        2        2
Eubacterium hallii et rel.              10      102       39       31       28
Lachnospira pectinoschiza et rel.       54      141      118       55       67
Burkholderia                             1        1        2        0        2
                                  Sample-6 Sample-7 Sample-8 Sample-9 Sample-10
Akkermansia                             16       26       30       19       125
Atopobium                                0        1        0        0         0
Megasphaera elsdenii et rel.            76       49       71      356       382
Oceanospirillum                          1        1        1        1         2
Prevotella tannerae et rel.             18       17       18        3       111
Catenibacterium mitsuokai et rel.        8        7        1        2         1
Helicobacter                             2        2        2        2         2
Eubacterium hallii et rel.              40       58       26       20        28
Lachnospira pectinoschiza et rel.       64       64       80       40        57
Burkholderia                             1        1        1        0         1
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
Akkermansia                        18       97       67
Atopobium                           0        1        0
Megasphaera elsdenii et rel.        3       14       46
Oceanospirillum                     1        2        1
Prevotella tannerae et rel.        63       19       16
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
Akkermansia                  0.0021279111 4.600209e-03 2.319704e-03
Atopobium                    0.0000000000 4.742483e-05 0.000000e+00
Megasphaera elsdenii et rel. 0.0003546519 6.639476e-04 1.592632e-03
Oceanospirillum              0.0001182173 9.484966e-05 3.462244e-05
Prevotella tannerae et rel.  0.0074476889 9.010718e-04 5.539591e-04
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
Akkermansia                  0.052595 0.077374 0.054917 0.159111
Atopobium                    0.000000 0.007848 0.000000 0.009903
Megasphaera elsdenii et rel. 0.021464 0.029370 0.045497 0.055163
Oceanospirillum              0.012391 0.011099 0.006706 0.014005
Prevotella tannerae et rel.  0.098511 0.034216 0.026826 0.038362
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
Akkermansia                  0.0021279111 4.600209e-03 2.319704e-03
Atopobium                    0.0000000000 4.742483e-05 0.000000e+00
Megasphaera elsdenii et rel. 0.0003546519 6.639476e-04 1.592632e-03
Oceanospirillum              0.0001182173 9.484966e-05 3.462244e-05
Prevotella tannerae et rel.  0.0074476889 9.010718e-04 5.539591e-04
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
                               Sample-1    Sample-2   Sample-3
Akkermansia                  -0.7333167  0.65311419  0.3442596
Atopobium                    -0.3664071  2.35344701 -0.3664071
Megasphaera elsdenii et rel. -1.0345716 -0.04063187  0.8182075
Oceanospirillum              -0.5570347  0.25304450 -0.5570347
Prevotella tannerae et rel.   0.6612837 -0.55293969 -0.7225946
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
                               Sample-1    Sample-2    Sample-3
Akkermansia                   0.4182158  0.92846743  0.72038722
Atopobium                    -1.1517915 -0.89241875 -1.18046753
Megasphaera elsdenii et rel. -0.4126041  0.05030439  0.55399343
Oceanospirillum              -0.7821978 -0.70271169 -0.86821022
Prevotella tannerae et rel.   1.0657708  0.18490370  0.09587261
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
Akkermansia                  1.278754 1.9912261 1.832509
Atopobium                    0.000000 0.3010300 0.000000
Megasphaera elsdenii et rel. 0.602060 1.1760913 1.672098
Oceanospirillum              0.301030 0.4771213 0.301030
Prevotella tannerae et rel.  1.806180 1.3010300 1.230449
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
Akkermansia                  1.278754 1.9912261 1.832509
Atopobium                    0.000000 0.3010300 0.000000
Megasphaera elsdenii et rel. 0.602060 1.1760913 1.672098
Oceanospirillum              0.301030 0.4771213 0.301030
Prevotella tannerae et rel.  1.806180 1.3010300 1.230449
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
Akkermansia                   1.2641691  2.2268092  1.7902914
Atopobium                    -3.5554744 -2.0404856 -3.1149834
Megasphaera elsdenii et rel. -0.4880348  0.3131385  1.4176161
Oceanospirillum              -1.4976474 -1.4908696 -2.0163711
Prevotella tannerae et rel.   2.5111514  0.6118107  0.3815242
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
Akkermansia                        19       98       68
Atopobium                           1        2        1
Megasphaera elsdenii et rel.        4       15       47
Oceanospirillum                     2        3        2
Prevotella tannerae et rel.        64       20       17
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
Akkermansia                        18       97       67
Atopobium                           0        1        0
Megasphaera elsdenii et rel.        3       14       46
Oceanospirillum                     1        2        1
Prevotella tannerae et rel.        63       19       16
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


