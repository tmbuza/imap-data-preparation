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
                                      Sample-1 Sample-2 Sample-3 Sample-4
Clostridium (sensu stricto)                 15       51       42      154
Clostridium symbiosum et rel.              220      247      214      111
Uncultured Clostridiales II                 27      205      284      155
Enterobacter aerogenes et rel.               7       10       14       36
Clostridium colinum et rel.                 24       30       11       10
Anaerovorax odorimutans et rel.              6       73       30       60
Phascolarctobacterium faecium et rel.        4        5       13        4
Eubacterium rectale et rel.                190       19       64       15
Peptostreptococcus micros et rel.            2        3        2        3
Dorea formicigenerans et rel.               42      267      419      135
                                      Sample-5 Sample-6 Sample-7 Sample-8
Clostridium (sensu stricto)                 21       39       42       21
Clostridium symbiosum et rel.              155      274      211      226
Uncultured Clostridiales II                 30       35      175       45
Enterobacter aerogenes et rel.              15       65        9        7
Clostridium colinum et rel.                 86        9        7       29
Anaerovorax odorimutans et rel.             12        9       45       14
Phascolarctobacterium faecium et rel.       33       20        8       11
Eubacterium rectale et rel.                 63      227       25      187
Peptostreptococcus micros et rel.            2        2        2        2
Dorea formicigenerans et rel.              189      331      480      156
                                      Sample-9 Sample-10
Clostridium (sensu stricto)                 24       104
Clostridium symbiosum et rel.              153       168
Uncultured Clostridiales II                 51       141
Enterobacter aerogenes et rel.              61        24
Clostridium colinum et rel.                  7        22
Anaerovorax odorimutans et rel.             14        32
Phascolarctobacterium faecium et rel.       21       103
Eubacterium rectale et rel.                 27        21
Peptostreptococcus micros et rel.            2         2
Dorea formicigenerans et rel.              101        70
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
Clostridium (sensu stricto)          15       51       42
Clostridium symbiosum et rel.       220      247      214
Uncultured Clostridiales II          27      205      284
Enterobacter aerogenes et rel.        7       10       14
Clostridium colinum et rel.          24       30       11
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
Clostridium (sensu stricto)    0.001773259 0.0024186664 0.0014541426
Clostridium symbiosum et rel.  0.026007802 0.0117139334 0.0074092026
Uncultured Clostridiales II    0.003191867 0.0097220905 0.0098327736
Enterobacter aerogenes et rel. 0.000827521 0.0004742483 0.0004847142
Clostridium colinum et rel.    0.002837215 0.0014227449 0.0003808469
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
Clostridium (sensu stricto)    0.048009 0.056077 0.043472 0.123198
Clostridium symbiosum et rel.  0.184841 0.123660 0.098255 0.104520
Uncultured Clostridiales II    0.064431 0.112608 0.113250 0.123600
Enterobacter aerogenes et rel. 0.032790 0.024821 0.025093 0.059450
Clostridium colinum et rel.    0.060741 0.043000 0.022242 0.031320
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
Clostridium (sensu stricto)    0.001773259 0.0024186664 0.0014541426
Clostridium symbiosum et rel.  0.026007802 0.0117139334 0.0074092026
Uncultured Clostridiales II    0.003191867 0.0097220905 0.0098327736
Enterobacter aerogenes et rel. 0.000827521 0.0004742483 0.0004847142
Clostridium colinum et rel.    0.002837215 0.0014227449 0.0003808469
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
Clostridium (sensu stricto)    -0.6651233  1.73874887  1.3511541
Clostridium symbiosum et rel.   0.2565546  0.40691215  0.2206504
Uncultured Clostridiales II    -0.5712673  1.65878223  2.0215188
Enterobacter aerogenes et rel. -0.4834495 -0.05597094  0.3603677
Clostridium colinum et rel.     0.3285117  0.56820398 -0.4893285
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
Clostridium (sensu stricto)     0.32658334  0.6319638  0.51392329
Clostridium symbiosum et rel.   1.72656831  1.3628714  1.23896230
Uncultured Clostridiales II     0.62497665  1.2760560  1.36593401
Enterobacter aerogenes et rel. -0.04301037 -0.0948094  0.03948761
Clostridium colinum et rel.     0.56454854  0.3899523 -0.06103679
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
Clostridium (sensu stricto)    1.204120 1.716003 1.633468
Clostridium symbiosum et rel.  2.344392 2.394452 2.332438
Uncultured Clostridiales II    1.447158 2.313867 2.454845
Enterobacter aerogenes et rel. 0.903090 1.041393 1.176091
Clostridium colinum et rel.    1.397940 1.491362 1.079181
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
Clostridium (sensu stricto)    1.204120 1.716003 1.633468
Clostridium symbiosum et rel.  2.344392 2.394452 2.332438
Uncultured Clostridiales II    1.447158 2.313867 2.454845
Enterobacter aerogenes et rel. 0.903090 1.041393 1.176091
Clostridium colinum et rel.    1.397940 1.491362 1.079181
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
Clostridium (sensu stricto)    1.0834602  1.58729959 1.32766789
Clostridium symbiosum et rel.  3.7599879  3.15920718 2.94647355
Uncultured Clostridiales II    1.6669407  2.97313111 3.22889707
Enterobacter aerogenes et rel. 0.3323087 -0.01322079 0.25231246
Clostridium colinum et rel.    1.5498317  1.06163354 0.02051085
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
Clostridium (sensu stricto)          16       52       43
Clostridium symbiosum et rel.       221      248      215
Uncultured Clostridiales II          28      206      285
Enterobacter aerogenes et rel.        8       11       15
Clostridium colinum et rel.          25       31       12
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
Clostridium (sensu stricto)          15       51       42
Clostridium symbiosum et rel.       220      247      214
Uncultured Clostridiales II          27      205      284
Enterobacter aerogenes et rel.        7       10       14
Clostridium colinum et rel.          24       30       11
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


