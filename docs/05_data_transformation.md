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
Actinomycetaceae                         0        1        0        1        0
Aerococcus                               0        0        0        0        0
Aeromonas                                0        0        0        0        0
Akkermansia                             18       97       67      256       21
Alcaligenes faecalis et rel.             1        2        3        2        2
Allistipes et rel.                     336       63       36       96       49
Anaerobiospirillum                       0        0        0        0        0
Anaerofustis                             0        1        0        0        0
Anaerostipes caccae et rel.            244      137       27       36       23
Anaerotruncus colihominis et rel.       12      108      203       68       15
                                  Sample-6 Sample-7 Sample-8 Sample-9 Sample-10
Actinomycetaceae                         0        0        0        0         0
Aerococcus                               0        0        0        0         0
Aeromonas                                0        0        0        0         0
Akkermansia                             16       26       30       19       125
Alcaligenes faecalis et rel.             2        2        2        2         7
Allistipes et rel.                      17       47       49       14        31
Anaerobiospirillum                       0        0        0        0         0
Anaerofustis                             0        0        0        0         0
Anaerostipes caccae et rel.             29       58       49       23        29
Anaerotruncus colihominis et rel.       36       31       40       19       720
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
Actinomycetaceae                    0        1        0
Aerococcus                          0        0        0
Aeromonas                           0        0        0
Akkermansia                        18       97       67
Alcaligenes faecalis et rel.        1        2        3
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
```

```r
otu_table(ps_rel)[1:5, 1:3]
```

```
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                                 Sample-1     Sample-2     Sample-3
Actinomycetaceae             0.0000000000 4.742483e-05 0.0000000000
Aerococcus                   0.0000000000 0.000000e+00 0.0000000000
Aeromonas                    0.0000000000 0.000000e+00 0.0000000000
Akkermansia                  0.0021279111 4.600209e-03 0.0023197036
Alcaligenes faecalis et rel. 0.0001182173 9.484966e-05 0.0001038673
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
Actinomycetaceae             0.000000 0.007848 0.000000 0.009903
Aerococcus                   0.000000 0.000000 0.000000 0.000000
Aeromonas                    0.000000 0.000000 0.000000 0.000000
Akkermansia                  0.052595 0.077374 0.054917 0.159111
Alcaligenes faecalis et rel. 0.012391 0.011099 0.011615 0.014005
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
Actinomycetaceae             0.0000000000 4.742483e-05 0.0000000000
Aerococcus                   0.0000000000 0.000000e+00 0.0000000000
Aeromonas                    0.0000000000 0.000000e+00 0.0000000000
Akkermansia                  0.0021279111 4.600209e-03 0.0023197036
Alcaligenes faecalis et rel. 0.0001182173 9.484966e-05 0.0001038673
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
Actinomycetaceae             -0.27806481  3.58008443 -0.27806481
Aerococcus                    0.00000000  0.00000000  0.00000000
Aeromonas                    -0.09513127 -0.09513127 -0.09513127
Akkermansia                  -0.73331670  0.65311419  0.34425962
Alcaligenes faecalis et rel. -0.93076516 -0.44189484 -0.09503583
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
Actinomycetaceae             -1.1517915 -0.8924188 -1.1804675
Aerococcus                   -1.1517915 -1.2167251 -1.1804675
Aeromonas                    -1.1517915 -1.2167251 -1.1804675
Akkermansia                   0.4182158  0.9284674  0.7203872
Alcaligenes faecalis et rel. -0.7821978 -0.7027117 -0.5559529
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
Actinomycetaceae             0.000000 0.3010300 0.000000
Aerococcus                   0.000000 0.0000000 0.000000
Aeromonas                    0.000000 0.0000000 0.000000
Akkermansia                  1.278754 1.9912261 1.832509
Alcaligenes faecalis et rel. 0.301030 0.4771213 0.602060
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
Actinomycetaceae             0.000000 0.3010300 0.000000
Aerococcus                   0.000000 0.0000000 0.000000
Aeromonas                    0.000000 0.0000000 0.000000
Akkermansia                  1.278754 1.9912261 1.832509
Alcaligenes faecalis et rel. 0.301030 0.4771213 0.602060
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
                              Sample-1  Sample-2  Sample-3
Actinomycetaceae             -3.555474 -2.040486 -3.114983
Aerococcus                   -3.555474 -3.359449 -3.114983
Aeromonas                    -3.555474 -3.359449 -3.114983
Akkermansia                   1.264169  2.226809  1.790291
Alcaligenes faecalis et rel. -1.497647 -1.490870 -1.169073
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
Actinomycetaceae                    1        2        1
Aerococcus                          1        1        1
Aeromonas                           1        1        1
Akkermansia                        19       98       68
Alcaligenes faecalis et rel.        2        3        4
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
Actinomycetaceae                    0        1        0
Aerococcus                          0        0        0
Aeromonas                           0        0        0
Akkermansia                        18       97       67
Alcaligenes faecalis et rel.        1        2        3
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


