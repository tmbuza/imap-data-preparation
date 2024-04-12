# (PART) DIMENSIONAL REDUCTION {-}

# Dimensionality Reduction {#data-reduction}
Data dimensional reduction is a crucial step in microbiome analysis, aimed at reducing the complexity of datasets while retaining important information. In this section, we'll explore techniques for reducing the number of features or variables in a phyloseq object, allowing for more efficient analysis and visualization of microbial community data.




## Prerequisite
- Preprocessed R objects
- Preprocessed phyloseq objects

> Either of the R or phyloseq objects can be used as input data.

## Import libraries and data objects

```r
set.seed(110912)

library(tidyverse, suppressPackageStartupMessages())
library(phyloseq)
library(microbiome)

cat("\nSaved phyloseq objects and dataframes\n\n")
```

```

Saved phyloseq objects and dataframes
```

```r
load("../imap-data-preparation/data/external/external_ps_objects.rda", verbose = T)
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



```r
library(phyloseq)
library(microbiome)

ps_raw <- ps_dietswap
ps_raw
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
```


