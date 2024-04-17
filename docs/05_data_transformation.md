# (PART) DATA TRANSFORMATION {-}

# Data Tranformation from Phyloseq Objects {#data-tarnsformation}

The process of data transformation in microbiome analysis involves converting raw or relative abundance values into matrices that are suitable for further analysis. Data transformation focuses on modifying the data itself. In this context, we will explore various methods used for this transformation on a phyloseq object.





## Import processed phyloseq object

```r
load("../imap-data-preparation/data/phyloseq_raw_rel_psextra_df_objects.rda", verbose = TRUE) 
```

```
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
```

```
OTU Table:          [10 taxa and 10 samples]
                     taxa are rows
                            Sample-1 Sample-2 Sample-3 Sample-4 Sample-5
Sporobacter termitidis            87     1187      713      764      138
Eubacterium hallii                10      102       39       31       28
Ruminococcus gnavus               15       56       65       41       47
Klebisiella pneumoniae             3        6        6       22        6
Veillonella                        8       19        6        5        7
Oceanospirillum                    1        2        1        2        1
Ruminococcus lactaris              3        9       23       39        7
Weissella                          1      180        2       13        2
Eggerthella lenta                  3       10        8       18        4
Uncultured Selenomonadaceae        0        0        0        0        0
                            Sample-6 Sample-7 Sample-8 Sample-9 Sample-10
Sporobacter termitidis            81      714      258      159       353
Eubacterium hallii                40       58       26       20        28
Ruminococcus gnavus               61      101       19       21        19
Klebisiella pneumoniae            32        4        5       51        20
Veillonella                        8       11       11       71        74
Oceanospirillum                    1        1        1        1         2
Ruminococcus lactaris             15       22        8        4         4
Weissella                          8        6        1       20         4
Eggerthella lenta                  7       13        6        6         8
Uncultured Selenomonadaceae        0        0        0        0         0
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
Sporobacter termitidis       87     1187      713
Eubacterium hallii           10      102       39
Ruminococcus gnavus          15       56       65
Klebisiella pneumoniae        3        6        6
Veillonella                   8       19        6
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
Sporobacter termitidis 0.0102849037 0.0562932752 0.0246858013
Eubacterium hallii     0.0011821728 0.0048373328 0.0013502752
Ruminococcus gnavus    0.0017732593 0.0026557906 0.0022504587
Klebisiella pneumoniae 0.0003546519 0.0002845490 0.0002077347
Veillonella            0.0009457383 0.0009010718 0.0002077347
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
Sporobacter termitidis 0.115836 0.273804 0.180029 0.277249
Eubacterium hallii     0.039194 0.079347 0.041890 0.055163
Ruminococcus gnavus    0.048009 0.058765 0.054090 0.063450
Klebisiella pneumoniae 0.021464 0.019225 0.016427 0.046464
Veillonella            0.035055 0.034216 0.016427 0.022144
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
Sporobacter termitidis 0.0102849037 0.0562932752 0.0246858013
Eubacterium hallii     0.0011821728 0.0048373328 0.0013502752
Ruminococcus gnavus    0.0017732593 0.0026557906 0.0022504587
Klebisiella pneumoniae 0.0003546519 0.0002845490 0.0002077347
Veillonella            0.0009457383 0.0009010718 0.0002077347
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
Sporobacter termitidis -0.7794277  1.90302575  1.37827866
Eubacterium hallii     -0.4477855  2.45470119  1.22737899
Ruminococcus gnavus    -0.4074006  1.29211010  1.48822307
Klebisiella pneumoniae -0.7431343 -0.01469595 -0.01469595
Veillonella             0.6448912  1.68310449  0.31813426
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
Sporobacter termitidis  1.23557381  2.0958437  1.7796631
Eubacterium hallii      0.12679267  0.9517496  0.4813434
Ruminococcus gnavus     0.32658334  0.6749183  0.7069387
Klebisiella pneumoniae -0.41260408 -0.3062820 -0.3038504
Veillonella             0.01979284  0.1849037 -0.3038504
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
                        Sample-1 Sample-2 Sample-3
Sporobacter termitidis 1.9444827 3.074816 2.853698
Eubacterium hallii     1.0413927 2.012837 1.602060
Ruminococcus gnavus    1.2041200 1.755875 1.819544
Klebisiella pneumoniae 0.6020600 0.845098 0.845098
Veillonella            0.9542425 1.301030 0.845098
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
                        Sample-1 Sample-2 Sample-3
Sporobacter termitidis 1.9444827 3.074816 2.853698
Eubacterium hallii     1.0413927 2.012837 1.602060
Ruminococcus gnavus    1.2041200 1.755875 1.819544
Klebisiella pneumoniae 0.6020600 0.845098 0.845098
Veillonella            0.9542425 1.301030 0.845098
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
                         Sample-1   Sample-2  Sample-3
Sporobacter termitidis  2.8332848  4.7278340  4.148346
Eubacterium hallii      0.6828174  2.2768872  1.254464
Ruminococcus gnavus     1.0834602  1.6801909  1.760214
Klebisiella pneumoniae -0.4880348 -0.5008398 -0.550034
Veillonella             0.4632754  0.6118107 -0.550034
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
Sporobacter termitidis       88     1188      714
Eubacterium hallii           11      103       40
Ruminococcus gnavus          16       57       66
Klebisiella pneumoniae        4        7        7
Veillonella                   9       20        7
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
Sporobacter termitidis       87     1187      713
Eubacterium hallii           10      102       39
Ruminococcus gnavus          15       56       65
Klebisiella pneumoniae        3        6        6
Veillonella                   8       19        6
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


