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
Ruminococcus obeum et rel.               89      901      620      476      457
Escherichia coli et rel.                  9       13       12       33       11
Ruminococcus callidus et rel.            17       78      124       38       40
Clostridium colinum et rel.              24       30       11       10       86
Anaerovorax odorimutans et rel.           6       73       30       60       12
Burkholderia                              1        1        2        0        2
Dorea formicigenerans et rel.            42      267      419      135      189
Coprobacillus catenaformis et rel.        3        2        3        4        4
Propionibacterium                         1        1        1        1        1
Granulicatella                            0        5        0        0        0
                                   Sample-6 Sample-7 Sample-8 Sample-9
Ruminococcus obeum et rel.              375      661      444      339
Escherichia coli et rel.                 54        8       17       73
Ruminococcus callidus et rel.            64       29       84       35
Clostridium colinum et rel.               9        7       29        7
Anaerovorax odorimutans et rel.           9       45       14       14
Burkholderia                              1        1        1        0
Dorea formicigenerans et rel.           331      480      156      101
Coprobacillus catenaformis et rel.        3        4        2        2
Propionibacterium                         1        1        1        1
Granulicatella                            1        0        0       16
                                   Sample-10
Ruminococcus obeum et rel.               599
Escherichia coli et rel.                  14
Ruminococcus callidus et rel.            545
Clostridium colinum et rel.               22
Anaerovorax odorimutans et rel.           32
Burkholderia                               1
Dorea formicigenerans et rel.             70
Coprobacillus catenaformis et rel.         2
Propionibacterium                          1
Granulicatella                             0
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
Ruminococcus obeum et rel.            89      901      620
Escherichia coli et rel.               9       13       12
Ruminococcus callidus et rel.         17       78      124
Clostridium colinum et rel.           24       30       11
Anaerovorax odorimutans et rel.        6       73       30
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
Ruminococcus obeum et rel.      0.0105213382 0.0427297733 0.0214659142
Escherichia coli et rel.        0.0010639556 0.0006165228 0.0004154693
Ruminococcus callidus et rel.   0.0020096938 0.0036991369 0.0042931828
Clostridium colinum et rel.     0.0028372148 0.0014227449 0.0003808469
Anaerovorax odorimutans et rel. 0.0007093037 0.0034620127 0.0010386733
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
Ruminococcus obeum et rel.      0.117166 0.237814 0.167759 0.217764
Escherichia coli et rel.        0.037182 0.028301 0.023232 0.056916
Ruminococcus callidus et rel.   0.051112 0.069370 0.074742 0.061081
Clostridium colinum et rel.     0.060741 0.043000 0.022242 0.031320
Anaerovorax odorimutans et rel. 0.030357 0.067106 0.036737 0.076780
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
Ruminococcus obeum et rel.      0.0105213382 0.0427297733 0.0214659142
Escherichia coli et rel.        0.0010639556 0.0006165228 0.0004154693
Ruminococcus callidus et rel.   0.0020096938 0.0036991369 0.0042931828
Clostridium colinum et rel.     0.0028372148 0.0014227449 0.0003808469
Anaerovorax odorimutans et rel. 0.0007093037 0.0034620127 0.0010386733
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
Ruminococcus obeum et rel.      -0.7310885  2.1112033  1.6508691
Escherichia coli et rel.        -0.6024728 -0.3248101 -0.3859653
Ruminococcus callidus et rel.   -0.8536852  0.8144821  1.3320113
Clostridium colinum et rel.      0.3285117  0.5682040 -0.4893285
Anaerovorax odorimutans et rel. -1.1540729  3.3528861  1.6899736
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
                                   Sample-1   Sample-2    Sample-3
Ruminococcus obeum et rel.       1.24755658 1.96698523  1.71679611
Escherichia coli et rel.         0.07597223 0.01802434 -0.02497819
Ruminococcus callidus et rel.    0.38938655 0.82763115  0.99464951
Clostridium colinum et rel.      0.56454854 0.38995231 -0.06103679
Anaerovorax odorimutans et rel. -0.11421078 0.79704017  0.36651647
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
Ruminococcus obeum et rel.      1.954243 2.955207 2.793092
Escherichia coli et rel.        1.000000 1.146128 1.113943
Ruminococcus callidus et rel.   1.255273 1.897627 2.096910
Clostridium colinum et rel.     1.397940 1.491362 1.079181
Anaerovorax odorimutans et rel. 0.845098 1.869232 1.491362
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
Ruminococcus obeum et rel.      1.954243 2.955207 2.793092
Escherichia coli et rel.        1.000000 1.146128 1.113943
Ruminococcus callidus et rel.   1.255273 1.897627 2.096910
Clostridium colinum et rel.     1.397940 1.491362 1.079181
Anaerovorax odorimutans et rel. 0.845098 1.869232 1.491362
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
                                 Sample-1  Sample-2   Sample-3
Ruminococcus obeum et rel.      2.8559753 4.4522524 4.00868942
Escherichia coli et rel.        0.5790592 0.2409833 0.10389246
Ruminococcus callidus et rel.   1.2074853 2.0097198 2.40246953
Clostridium colinum et rel.     1.5498317 1.0616335 0.02051085
Anaerovorax odorimutans et rel. 0.1815673 1.9437894 0.99589050
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
Ruminococcus obeum et rel.            90      902      621
Escherichia coli et rel.              10       14       13
Ruminococcus callidus et rel.         18       79      125
Clostridium colinum et rel.           25       31       12
Anaerovorax odorimutans et rel.        7       74       31
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
Ruminococcus obeum et rel.            89      901      620
Escherichia coli et rel.               9       13       12
Ruminococcus callidus et rel.         17       78      124
Clostridium colinum et rel.           24       30       11
Anaerovorax odorimutans et rel.        6       73       30
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


