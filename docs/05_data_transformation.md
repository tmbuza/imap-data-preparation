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
                           Sample-1 Sample-2 Sample-3 Sample-4 Sample-5
Parabacteroides distasonis      187       28       17       19       38
Eggerthella lenta                 3       10        8       18        4
Coprococcus eutactus             49      355      477      329      211
Serratia                          1        1        0        3        0
Clostridium felsineum             0        0        0        0        0
Clostridium orbiscindens         85      750      510      691      265
Anaerovorax odorimutans           6       73       30       60       12
Bifidobacterium                  43       25      183      493       22
Fusobacteria                      4        5        5        5        5
Eubacterium limosum               1        1        1        1        1
                           Sample-6 Sample-7 Sample-8 Sample-9 Sample-10
Parabacteroides distasonis       12       27       19       10        15
Eggerthella lenta                 7       13        6        6         8
Coprococcus eutactus            137      332      267      146        52
Serratia                          9        0        1       16         2
Clostridium felsineum             0        0        0        0         0
Clostridium orbiscindens        229      671      243      447       149
Anaerovorax odorimutans           9       45       14       14        32
Bifidobacterium                  24      116       62      151        20
Fusobacteria                      8        5        5        4         5
Eubacterium limosum               1        1        1        1         1
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
Parabacteroides distasonis      187       28       17
Eggerthella lenta                 3       10        8
Coprococcus eutactus             49      355      477
Serratia                          1        1        0
Clostridium felsineum             0        0        0
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
Parabacteroides distasonis 0.0221066320 1.327895e-03 0.0005885815
Eggerthella lenta          0.0003546519 4.742483e-04 0.0002769795
Coprococcus eutactus       0.0057926469 1.683582e-02 0.0165149050
Serratia                   0.0001182173 4.742483e-05 0.0000000000
Clostridium felsineum      0.0000000000 0.000000e+00 0.0000000000
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
Parabacteroides distasonis 0.170268 0.041541 0.027652 0.043177
Eggerthella lenta          0.021464 0.024821 0.018968 0.042025
Coprococcus eutactus       0.086847 0.148417 0.146985 0.180595
Serratia                   0.012391 0.007848 0.000000 0.017152
Clostridium felsineum      0.000000 0.000000 0.000000 0.000000
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
Parabacteroides distasonis 0.0221066320 1.327895e-03 0.0005885815
Eggerthella lenta          0.0003546519 4.742483e-04 0.0002769795
Coprococcus eutactus       0.0057926469 1.683582e-02 0.0165149050
Serratia                   0.0001182173 4.742483e-05 0.0000000000
Clostridium felsineum      0.0000000000 0.000000e+00 0.0000000000
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
Parabacteroides distasonis  0.994190116 -0.7348082 -1.1759727
Eggerthella lenta          -0.821628886  2.3094335  1.6883264
Coprococcus eutactus       -0.006970883  2.0022117  2.3038387
Serratia                    0.033901004  0.0339010 -0.6285507
Clostridium felsineum       0.000000000  0.0000000  0.0000000
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
Parabacteroides distasonis  1.6403370  0.3587491  0.1216220
Eggerthella lenta          -0.4126041 -0.0948094 -0.1906353
Coprococcus eutactus        0.9341423  1.5320092  1.5988932
Serratia                   -0.7821978 -0.8924188 -1.1804675
Clostridium felsineum      -1.1517915 -1.2167251 -1.1804675
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
                           Sample-1 Sample-2  Sample-3
Parabacteroides distasonis 2.274158 1.462398 1.2552725
Eggerthella lenta          0.602060 1.041393 0.9542425
Coprococcus eutactus       1.698970 2.551450 2.6794279
Serratia                   0.301030 0.301030 0.0000000
Clostridium felsineum      0.000000 0.000000 0.0000000
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
                           Sample-1 Sample-2  Sample-3
Parabacteroides distasonis 2.274158 1.462398 1.2552725
Eggerthella lenta          0.602060 1.041393 0.9542425
Coprococcus eutactus       1.698970 2.551450 2.6794279
Serratia                   0.301030 0.301030 0.0000000
Clostridium felsineum      0.000000 0.000000 0.0000000
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
                             Sample-1    Sample-2   Sample-3
Parabacteroides distasonis  3.5975863  0.99349896  0.4403647
Eggerthella lenta          -0.4880348 -0.01322079 -0.2817700
Coprococcus eutactus        2.2604993  3.52148760  3.7467280
Serratia                   -1.4976474 -2.04048558 -3.1149834
Clostridium felsineum      -3.5554744 -3.35944895 -3.1149834
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
Parabacteroides distasonis      188       29       18
Eggerthella lenta                 4       11        9
Coprococcus eutactus             50      356      478
Serratia                          2        2        1
Clostridium felsineum             1        1        1
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
Parabacteroides distasonis      187       28       17
Eggerthella lenta                 3       10        8
Coprococcus eutactus             49      355      477
Serratia                          1        1        0
Clostridium felsineum             0        0        0
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


