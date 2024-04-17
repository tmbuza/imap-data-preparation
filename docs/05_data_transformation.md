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
Ruminococcus bromii et rel.            72      403      106       60       30
Enterococcus                            3      232        5        5        4
Lactococcus                             1        4        2       81        1
Roseburia intestinalis et rel.         34       10       18        7       24
Bacteroides fragilis et rel.          443       21       73       29       33
Bacteroides intestinalis et rel.       12        1        3        6        3
Prevotella oralis et rel.              22      747     1832       10     1505
Eubacterium ventriosum et rel.         10       20       13       11       14
Uncultured Clostridiales II            27      205      284      155       30
Granulicatella                          0        5        0        0        0
                                 Sample-6 Sample-7 Sample-8 Sample-9 Sample-10
Ruminococcus bromii et rel.            19       17      189       16        55
Enterococcus                           11        4        6      692         6
Lactococcus                             6        1        1        7         1
Roseburia intestinalis et rel.         31       14       70       13         8
Bacteroides fragilis et rel.           14       22      560        8        26
Bacteroides intestinalis et rel.        1        1        2        0         1
Prevotella oralis et rel.            1454       75     1806      138       637
Eubacterium ventriosum et rel.         18       19       25       10        15
Uncultured Clostridiales II            35      175       45       51       141
Granulicatella                          1        0        0       16         0
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
Ruminococcus bromii et rel.          72      403      106
Enterococcus                          3      232        5
Lactococcus                           1        4        2
Roseburia intestinalis et rel.       34       10       18
Bacteroides fragilis et rel.        443       21       73
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
Ruminococcus bromii et rel.    0.0085116444 0.0191122072 3.669979e-03
Enterococcus                   0.0003546519 0.0110025609 1.731122e-04
Lactococcus                    0.0001182173 0.0001896993 6.924488e-05
Roseburia intestinalis et rel. 0.0040193876 0.0004742483 6.232040e-04
Bacteroides fragilis et rel.   0.0523702565 0.0009959215 2.527438e-03
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
Ruminococcus bromii et rel.    0.105337 0.158212 0.069095 0.076780
Enterococcus                   0.021464 0.119828 0.014995 0.022144
Lactococcus                    0.012391 0.015697 0.009484 0.089241
Roseburia intestinalis et rel. 0.072315 0.024821 0.028454 0.026203
Bacteroides fragilis et rel.   0.263854 0.035973 0.057326 0.053352
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
Ruminococcus bromii et rel.    0.0085116444 0.0191122072 3.669979e-03
Enterococcus                   0.0003546519 0.0110025609 1.731122e-04
Lactococcus                    0.0001182173 0.0001896993 6.924488e-05
Roseburia intestinalis et rel. 0.0040193876 0.0004742483 6.232040e-04
Bacteroides fragilis et rel.   0.0523702565 0.0009959215 2.527438e-03
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
                                 Sample-1   Sample-2    Sample-3
Ruminococcus bromii et rel.     0.8452788  2.4252435  1.19837397
Enterococcus                   -0.5423550  5.8629260  0.09658265
Lactococcus                    -0.2903369  2.3306971  0.86948912
Roseburia intestinalis et rel.  1.0432030 -0.4633681  0.24802758
Bacteroides fragilis et rel.    1.3124945 -1.1679540 -0.16660347
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
Ruminococcus bromii et rel.     1.1359287  1.5911880  0.9246048
Enterococcus                   -0.4126041  1.3336806 -0.3732941
Lactococcus                    -0.7821978 -0.4637090 -0.6855514
Roseburia intestinalis et rel.  0.7439592 -0.0948094  0.1459789
Bacteroides fragilis et rel.    2.0985693  0.2294970  0.7584797
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
                               Sample-1 Sample-2  Sample-3
Ruminococcus bromii et rel.    1.863323 2.606381 2.0293838
Enterococcus                   0.602060 2.367356 0.7781513
Lactococcus                    0.301030 0.698970 0.4771213
Roseburia intestinalis et rel. 1.544068 1.041393 1.2787536
Bacteroides fragilis et rel.   2.647383 1.342423 1.8692317
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
                               Sample-1 Sample-2  Sample-3
Ruminococcus bromii et rel.    1.863323 2.606381 2.0293838
Enterococcus                   0.602060 2.367356 0.7781513
Lactococcus                    0.301030 0.698970 0.4771213
Roseburia intestinalis et rel. 1.544068 1.041393 1.2787536
Bacteroides fragilis et rel.   2.647383 1.342423 1.8692317
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
                                 Sample-1    Sample-2   Sample-3
Ruminococcus bromii et rel.     2.6443928  3.64818403  2.2463088
Enterococcus                   -0.4880348  3.09665162 -0.7170881
Lactococcus                    -1.4976474 -0.87803411 -1.5055455
Roseburia intestinalis et rel.  1.8963532 -0.01322079  0.4959345
Bacteroides fragilis et rel.    4.4595952  0.71009732  1.8754492
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
Ruminococcus bromii et rel.          73      404      107
Enterococcus                          4      233        6
Lactococcus                           2        5        3
Roseburia intestinalis et rel.       35       11       19
Bacteroides fragilis et rel.        444       22       74
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
Ruminococcus bromii et rel.          72      403      106
Enterococcus                          3      232        5
Lactococcus                           1        4        2
Roseburia intestinalis et rel.       34       10       18
Bacteroides fragilis et rel.        443       21       73
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


