# (PART) DATA TRANSFORMATION {-}

# Data Tranformation from Phyloseq Objects {#data-tarnsformation}

The process of data transformation in microbiome analysis involves converting raw or relative abundance values into matrices that are suitable for further analysis. Data transformation focuses on modifying the data itself. In this context, we will explore various methods used for this transformation on a phyloseq object.





## Import processed phyloseq object

```r
load("../imap-data-preparation/data/phyloseq_raw_rel_psextra_df_objects.rda", verbose = TRUE) 
Loading objects:
  ps_raw
  ps_rel
  psextra_raw
  psextra_rel
  ps_df
```

## Raw Abundance


```r
ps_raw <- psextra_raw
otu_table(ps_raw)[1:10, 1:10]
OTU Table:          [10 taxa and 10 samples]
                     taxa are rows
                       Sample-1 Sample-2 Sample-3 Sample-4 Sample-5 Sample-6
Eubacterium limosum           1        1        1        1        1        1
Staphylococcus                0        0        0        0        0        0
Oceanospirillum               1        2        1        2        1        1
Ruminococcus obeum           89      901      620      476      457      375
Burkholderia                  1        1        2        0        2        1
Clostridium sphenoides       41      127      133       63       94       90
Fusobacteria                  4        5        5        5        5        8
Bifidobacterium              43       25      183      493       22       24
Methylobacterium              0        0        0        0        0        0
Clostridium ramosum           2        2        2        2        2        2
                       Sample-7 Sample-8 Sample-9 Sample-10
Eubacterium limosum           1        1        1         1
Staphylococcus                0        0        0         0
Oceanospirillum               1        1        1         2
Ruminococcus obeum          661      444      339       599
Burkholderia                  1        1        0         1
Clostridium sphenoides      153       92       55        84
Fusobacteria                  5        5        4         5
Bifidobacterium             116       62      151        20
Methylobacterium              0        0        0         0
Clostridium ramosum           3        2        2         2
```

##  No Transformation

> Similar to raw abundance


```r
(ps_identity <- microbiome::transform(ps_raw, 'identity'))
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]

cat("\n\n")
otu_table(ps_identity)[1:5, 1:3]
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                    Sample-1 Sample-2 Sample-3
Eubacterium limosum        1        1        1
Staphylococcus             0        0        0
Oceanospirillum            1        2        1
Ruminococcus obeum        89      901      620
Burkholderia               1        1        2
```

## Relative abundance


```r
(ps_rel = phyloseq::transform_sample_counts(ps_raw, function(x){x / sum(x)}))
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]

otu_table(ps_rel)[1:5, 1:3]
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                        Sample-1     Sample-2     Sample-3
Eubacterium limosum 0.0001182173 4.742483e-05 3.462244e-05
Staphylococcus      0.0000000000 0.000000e+00 0.000000e+00
Oceanospirillum     0.0001182173 9.484966e-05 3.462244e-05
Ruminococcus obeum  0.0105213382 4.272977e-02 2.146591e-02
Burkholderia        0.0001182173 4.742483e-05 6.924488e-05
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
OTU Table:          [5 taxa and 4 samples]
                     taxa are rows
                    Sample-1 Sample-2 Sample-3 Sample-4
Eubacterium limosum 0.012391 0.007848 0.006706 0.009903
Staphylococcus      0.000000 0.000000 0.000000 0.000000
Oceanospirillum     0.012391 0.011099 0.006706 0.014005
Ruminococcus obeum  0.117166 0.237814 0.167759 0.217764
Burkholderia        0.012391 0.007848 0.009484 0.000000
```


## Compositional Version
Compositional data represents relative proportions or percentages of different microbial taxa within a sample, rather than the absolute abundance of each taxon. This transformation is necessary because raw abundance data is typically affected by various factors, such as sequencing depth, which can lead to spurious correlations and biases in downstream analyses.


```r
(ps_compositional <- microbiome::transform(ps_raw, 'compositional'))
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]

cat("\n\n")
otu_table(ps_compositional)[1:5, 1:3]
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                        Sample-1     Sample-2     Sample-3
Eubacterium limosum 0.0001182173 4.742483e-05 3.462244e-05
Staphylococcus      0.0000000000 0.000000e+00 0.000000e+00
Oceanospirillum     0.0001182173 9.484966e-05 3.462244e-05
Ruminococcus obeum  0.0105213382 4.272977e-02 2.146591e-02
Burkholderia        0.0001182173 4.742483e-05 6.924488e-05
```

## Z-transform for OTUs


```r
(ps_z_otu <- microbiome::transform(ps_raw, 'Z', 'OTU'))
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]

cat("\n\n")
otu_table(ps_z_otu)[1:5, 1:3]
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                      Sample-1   Sample-2   Sample-3
Eubacterium limosum -0.2855358 -0.2855358 -0.2855358
Staphylococcus      -0.2685840 -0.2685840 -0.2685840
Oceanospirillum     -0.5570347  0.2530445 -0.5570347
Ruminococcus obeum  -0.7310885  2.1112033  1.6508691
Burkholderia         0.3632341  0.3632341  1.1488042
```

## Z-transform for Samples


```r
(ps_z_sample <- microbiome::transform(ps_raw, 'Z', 'sample'))
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]

cat("\n\n")
otu_table(ps_z_sample)[1:5, 1:3]
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                      Sample-1   Sample-2   Sample-3
Eubacterium limosum -0.7821978 -0.8924188 -0.8682102
Staphylococcus      -1.1517915 -1.2167251 -1.1804675
Oceanospirillum     -0.7821978 -0.7027117 -0.8682102
Ruminococcus obeum   1.2475566  1.9669852  1.7167961
Burkholderia        -0.7821978 -0.8924188 -0.6855514
```

## Log10 Transform


```r
(ps_log10 <- microbiome::transform(ps_raw, 'log10'))
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]

cat("\n\n")
otu_table(ps_log10)[1:5, 1:3]
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                    Sample-1  Sample-2  Sample-3
Eubacterium limosum 0.301030 0.3010300 0.3010300
Staphylococcus      0.000000 0.0000000 0.0000000
Oceanospirillum     0.301030 0.4771213 0.3010300
Ruminococcus obeum  1.954243 2.9552065 2.7930916
Burkholderia        0.301030 0.3010300 0.4771213
```

## Log10p Transform


```r
(ps_log10p <- microbiome::transform(ps_raw, 'log10p'))
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]

cat("\n\n")
otu_table(ps_log10p)[1:5, 1:3]
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                    Sample-1  Sample-2  Sample-3
Eubacterium limosum 0.301030 0.3010300 0.3010300
Staphylococcus      0.000000 0.0000000 0.0000000
Oceanospirillum     0.301030 0.4771213 0.3010300
Ruminococcus obeum  1.954243 2.9552065 2.7930916
Burkholderia        0.301030 0.3010300 0.4771213
```

## CLR Transform

- Note that small pseudocount is added if data contains zeroes

```r
(ps_clr <- microbiome::transform(ps_raw, 'clr'))
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]

cat("\n\n")
otu_table(ps_clr)[1:5, 1:3]
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                     Sample-1  Sample-2  Sample-3
Eubacterium limosum -1.497647 -2.040486 -2.016371
Staphylococcus      -3.555474 -3.359449 -3.114983
Oceanospirillum     -1.497647 -1.490870 -2.016371
Ruminococcus obeum   2.855975  4.452252  4.008689
Burkholderia        -1.497647 -2.040486 -1.505545
```


## Shift the baseline


```r
(ps_shift <- microbiome::transform(ps_raw, 'shift', shift=1))
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]

cat("\n\n")
otu_table(ps_shift)[1:5, 1:3]
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                    Sample-1 Sample-2 Sample-3
Eubacterium limosum        2        2        2
Staphylococcus             1        1        1
Oceanospirillum            2        3        2
Ruminococcus obeum        90      902      621
Burkholderia               2        2        3
```

## Data Scaling


```r
(ps_scale <- microbiome::transform(ps_raw, 'scale', scale=1))
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]

cat("\n\n")
otu_table(ps_scale)[1:5, 1:3]
OTU Table:          [5 taxa and 3 samples]
                     taxa are rows
                    Sample-1 Sample-2 Sample-3
Eubacterium limosum        1        1        1
Staphylococcus             0        0        0
Oceanospirillum            1        2        1
Ruminococcus obeum        89      901      620
Burkholderia               1        1        2
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


