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
Bacteroides fragilis et rel.          443       21       73       29       33
Actinomycetaceae                        0        1        0        1        0
Ruminococcus lactaris et rel.           3        9       23       39        7
Bryantella formatexigens et rel.       78      106      126       47       54
Uncultured Clostridiales I             19      120     1183      399       19
Eubacterium hallii et rel.             10      102       39       31       28
Ruminococcus obeum et rel.             89      901      620      476      457
Clostridium difficile et rel.          10       21       94      546       25
Enterobacter aerogenes et rel.          7       10       14       36       15
Asteroleplasma et rel.                  0        0        0        0        0
                                 Sample-6 Sample-7 Sample-8 Sample-9 Sample-10
Bacteroides fragilis et rel.           14       22      560        8        26
Actinomycetaceae                        0        0        0        0         0
Ruminococcus lactaris et rel.          15       22        8        4         4
Bryantella formatexigens et rel.       52       73       43       41        48
Uncultured Clostridiales I             43      137       70       59       315
Eubacterium hallii et rel.             40       58       26       20        28
Ruminococcus obeum et rel.            375      661      444      339       599
Clostridium difficile et rel.          45       45       41       65        82
Enterobacter aerogenes et rel.         65        9        7       61        24
Asteroleplasma et rel.                  0        0        0        0         0
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
Bacteroides fragilis et rel.          443       21       73
Actinomycetaceae                        0        1        0
Ruminococcus lactaris et rel.           3        9       23
Bryantella formatexigens et rel.       78      106      126
Uncultured Clostridiales I             19      120     1183
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
Bacteroides fragilis et rel.     0.0523702565 9.959215e-04 0.0025274383
Actinomycetaceae                 0.0000000000 4.742483e-05 0.0000000000
Ruminococcus lactaris et rel.    0.0003546519 4.268235e-04 0.0007963162
Bryantella formatexigens et rel. 0.0092209481 5.027032e-03 0.0043624277
Uncultured Clostridiales I       0.0022461284 5.690980e-03 0.0409583492
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
Bacteroides fragilis et rel.     0.263854 0.035973 0.057326 0.053352
Actinomycetaceae                 0.000000 0.007848 0.000000 0.009903
Ruminococcus lactaris et rel.    0.021464 0.023547 0.032165 0.061881
Bryantella formatexigens et rel. 0.109655 0.080891 0.075344 0.067940
Uncultured Clostridiales I       0.054038 0.086080 0.232740 0.199116
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
Bacteroides fragilis et rel.     0.0523702565 9.959215e-04 0.0025274383
Actinomycetaceae                 0.0000000000 4.742483e-05 0.0000000000
Ruminococcus lactaris et rel.    0.0003546519 4.268235e-04 0.0007963162
Bryantella formatexigens et rel. 0.0092209481 5.027032e-03 0.0043624277
Uncultured Clostridiales I       0.0022461284 5.690980e-03 0.0409583492
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
Bacteroides fragilis et rel.      1.3124945 -1.1679540 -0.1666035
Actinomycetaceae                 -0.2780648  3.5800844 -0.2780648
Ruminococcus lactaris et rel.    -0.6823416  0.4788310  1.5882717
Bryantella formatexigens et rel.  0.4035063  0.7474051  0.9416489
Uncultured Clostridiales I       -0.7313669  1.0472965  3.3010503
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
Bacteroides fragilis et rel.      2.0985693  0.2294970  0.7584797
Actinomycetaceae                 -1.1517915 -0.8924188 -1.1804675
Ruminococcus lactaris et rel.    -0.4126041 -0.1394027  0.2512205
Bryantella formatexigens et rel.  1.1780462  0.9695756  1.0018003
Uncultured Clostridiales I        0.4455659  1.0271063  2.0075089
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
                                 Sample-1 Sample-2 Sample-3
Bacteroides fragilis et rel.     2.647383 1.342423 1.869232
Actinomycetaceae                 0.000000 0.301030 0.000000
Ruminococcus lactaris et rel.    0.602060 1.000000 1.380211
Bryantella formatexigens et rel. 1.897627 2.029384 2.103804
Uncultured Clostridiales I       1.301030 2.082785 3.073352
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
                                 Sample-1 Sample-2 Sample-3
Bacteroides fragilis et rel.     2.647383 1.342423 1.869232
Actinomycetaceae                 0.000000 0.301030 0.000000
Ruminococcus lactaris et rel.    0.602060 1.000000 1.380211
Bryantella formatexigens et rel. 1.897627 2.029384 2.103804
Uncultured Clostridiales I       1.301030 2.082785 3.073352
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
Bacteroides fragilis et rel.      4.4595952  0.7100973  1.8754492
Actinomycetaceae                 -3.5554744 -2.0404856 -3.1149834
Ruminococcus lactaris et rel.    -0.4880348 -0.1146760  0.7351642
Bryantella formatexigens et rel.  2.7242794  2.3152189  2.4184061
Uncultured Clostridiales I        1.3178115  2.4388711  4.6543952
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
Bacteroides fragilis et rel.          444       22       74
Actinomycetaceae                        1        2        1
Ruminococcus lactaris et rel.           4       10       24
Bryantella formatexigens et rel.       79      107      127
Uncultured Clostridiales I             20      121     1184
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
Bacteroides fragilis et rel.          443       21       73
Actinomycetaceae                        0        1        0
Ruminococcus lactaris et rel.           3        9       23
Bryantella formatexigens et rel.       78      106      126
Uncultured Clostridiales I             19      120     1183
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


