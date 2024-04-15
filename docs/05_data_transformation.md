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
                                 Sample-1 Sample-2 Sample-3 Sample-4 Sample-5
Clostridium cellulosi et rel.          91      917      495      310      728
Oceanospirillum                         1        2        1        2        1
Bifidobacterium                        43       25      183      493       22
Clostridium sphenoides et rel.         41      127      133       63       94
Micrococcaceae                          0        0        0        0        0
Peptococcus niger et rel.               2        4       11       10        2
Burkholderia                            1        1        2        0        2
Streptococcus bovis et rel.           112       41       95       27      110
Lactobacillus salivarius et rel.        1        1        3       16        1
Mitsuokella multiacida et rel.          1        2       44        5        1
                                 Sample-6 Sample-7 Sample-8 Sample-9 Sample-10
Clostridium cellulosi et rel.         281     2528      992      491      1899
Oceanospirillum                         1        1        1        1         2
Bifidobacterium                        24      116       62      151        20
Clostridium sphenoides et rel.         90      153       92       55        84
Micrococcaceae                          0        0        0        0         0
Peptococcus niger et rel.               5        8       13        2         7
Burkholderia                            1        1        1        0         1
Streptococcus bovis et rel.           364       17       42      925        90
Lactobacillus salivarius et rel.        1        1        2        2         1
Mitsuokella multiacida et rel.        102       19       21      120       392
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
Clostridium cellulosi et rel.        91      917      495
Oceanospirillum                       1        2        1
Bifidobacterium                      43       25      183
Clostridium sphenoides et rel.       41      127      133
Micrococcaceae                        0        0        0
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
Clostridium cellulosi et rel.  0.0107577728 4.348857e-02 1.713811e-02
Oceanospirillum                0.0001182173 9.484966e-05 3.462244e-05
Bifidobacterium                0.0050833432 1.185621e-03 6.335907e-03
Clostridium sphenoides et rel. 0.0048469086 6.022954e-03 4.604785e-03
Micrococcaceae                 0.0000000000 0.000000e+00 0.000000e+00
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
Clostridium cellulosi et rel.  0.118481 0.239958 0.149754 0.175247
Oceanospirillum                0.012391 0.011099 0.006706 0.014005
Bifidobacterium                0.081344 0.039251 0.090839 0.221683
Clostridium sphenoides et rel. 0.079426 0.088561 0.077412 0.078680
Micrococcaceae                 0.000000 0.000000 0.000000 0.000000
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
Clostridium cellulosi et rel.  0.0107577728 4.348857e-02 1.713811e-02
Oceanospirillum                0.0001182173 9.484966e-05 3.462244e-05
Bifidobacterium                0.0050833432 1.185621e-03 6.335907e-03
Clostridium sphenoides et rel. 0.0048469086 6.022954e-03 4.604785e-03
Micrococcaceae                 0.0000000000 0.000000e+00 0.000000e+00
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
Clostridium cellulosi et rel.  -0.7468356  0.7886759  0.3777517
Oceanospirillum                -0.5570347  0.2530445 -0.5570347
Bifidobacterium                 0.1091362 -0.4247936  1.5611951
Clostridium sphenoides et rel. -0.5528961  1.0892951  1.1568029
Micrococcaceae                  0.0000000  0.0000000  0.0000000
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
Clostridium cellulosi et rel.   1.2592760  1.9752118  1.6155457
Oceanospirillum                -0.7821978 -0.7027117 -0.8682102
Bifidobacterium                 0.8659801  0.3076574  1.1688197
Clostridium sphenoides et rel.  0.8411751  1.0534195  1.0259704
Micrococcaceae                 -1.1517915 -1.2167251 -1.1804675
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
Clostridium cellulosi et rel.  1.963788 2.9628427 2.695482
Oceanospirillum                0.301030 0.4771213 0.301030
Bifidobacterium                1.643453 1.4149733 2.264818
Clostridium sphenoides et rel. 1.623249 2.1072100 2.127105
Micrococcaceae                 0.000000 0.0000000 0.000000
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
Clostridium cellulosi et rel.  1.963788 2.9628427 2.695482
Oceanospirillum                0.301030 0.4771213 0.301030
Bifidobacterium                1.643453 1.4149733 2.264818
Clostridium sphenoides et rel. 1.623249 2.1072100 2.127105
Micrococcaceae                 0.000000 0.0000000 0.000000
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
Clostridium cellulosi et rel.   2.878162  4.4698476  3.783731
Oceanospirillum                -1.497647 -1.4908696 -2.016371
Bifidobacterium                 2.130295  0.8817133  2.790378
Clostridium sphenoides et rel.  2.082832  2.4953993  2.472265
Micrococcaceae                 -3.555474 -3.3594489 -3.114983
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
Clostridium cellulosi et rel.        92      918      496
Oceanospirillum                       2        3        2
Bifidobacterium                      44       26      184
Clostridium sphenoides et rel.       42      128      134
Micrococcaceae                        1        1        1
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
Clostridium cellulosi et rel.        91      917      495
Oceanospirillum                       1        2        1
Bifidobacterium                      43       25      183
Clostridium sphenoides et rel.       41      127      133
Micrococcaceae                        0        0        0
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


