# (PART) QIIME2 OUTPUT {-}
# Prepare Qiime2 Output

## Data directories

```r
if (!dir.exists('data')) {dir.create('data')}
if (!dir.exists('data/qiime2')) {dir.create('data/qiime2')}

library(tidyverse, suppressPackageStartupMessages())
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.1.4     ✔ readr     2.1.5
## ✔ forcats   1.0.0     ✔ stringr   1.5.1
## ✔ ggplot2   3.4.4     ✔ tibble    3.2.1
## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
## ✔ purrr     1.0.2     
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

## Qiime2 metadata

```r
read_tsv("../imap-qiime2-bioinformatics/resources/metadata/qiime2_sample_metadata.tsv", show_col_types = FALSE) %>% 
  dplyr::rename(sample_id="sample-id") %>% 
  write_csv("data/qiime2/qiime2_tidy_metadata.csv")
```


## Qiime2 otutable

```r
read_tsv("../imap-qiime2-bioinformatics/qiime2_process/export/feature-table.tsv", skip = 1, show_col_types = FALSE) %>%
  dplyr::rename(feature='#OTU ID') %>%
  select(-starts_with('Mock')) %>% 
  mutate_at(2:ncol(.), as.numeric) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  pivot_longer(-feature, names_to = "sample_id", values_to = "count") %>% 
  relocate(sample_id, .before = feature) %>% 
  write_csv("data/qiime2/qiime2_tidy_otutable.csv")
```


## Qiime2 taxonomy

```r
read_tsv("../imap-qiime2-bioinformatics/qiime2_process/export/taxonomy.tsv", show_col_types=FALSE) %>% 
  dplyr::rename(feature="Feature ID") %>% 
  distinct() %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  mutate(Taxon = str_replace_all(Taxon, "; s__$", ""),
         Taxon = str_replace_all(Taxon, "; g__$", ""),
         Taxon = str_replace_all(Taxon, "; f__$", ""),
         Taxon = str_replace_all(Taxon, "; o__$", ""),
         Taxon = str_replace_all(Taxon, "; c__$", ""),
         Taxon = str_replace_all(Taxon, "; p__$", ""),
         Taxon = str_replace_all(Taxon, "; k__$", ""),
         Taxon = str_replace_all(Taxon, "\\[|\\]", ""),
         Taxon = str_replace_all(Taxon, "\\s", "")) %>%
  dplyr::filter(!grepl("s__*", Taxon)) %>%
  dplyr::filter(grepl("g__*", Taxon)) %>% 
  select(-Confidence) %>% 
  mutate(Taxon = str_replace_all(Taxon, "\\w__", "")) %>% 
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";") %>% 
  write_csv("data/qiime2/qiime2_tidy_taxonomy.csv")
```


## Qiime2 composite

```r
# QIIME2 composite
library(tidyverse, suppressPackageStartupMessages())

qiime2_tidy_metadata <- read_csv("data/qiime2/qiime2_tidy_metadata.csv", show_col_types = FALSE)
qiime2_tidy_otutable <- read_csv("data/qiime2/qiime2_tidy_otutable.csv", show_col_types = FALSE)
qiime2_tidy_taxonomy <- read_csv("data/qiime2/qiime2_tidy_taxonomy.csv", show_col_types = FALSE)

qiime2_composite <- inner_join(qiime2_tidy_metadata, qiime2_tidy_otutable, by = "sample_id") %>% 
  inner_join(., qiime2_tidy_taxonomy, by = "feature") %>% 
  group_by(sample_id) %>% 
  mutate(rel_abund = count/sum(count)) %>% 
  ungroup() %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  relocate(count, .before = rel_abund) 

write_csv(qiime2_composite, "data/qiime2/qiime2_composite.csv")
```

## Create a phyloseq object

```r
library(tidyverse)
library(phyloseq)
```

```
## Warning in .recacheSubclasses(def@className, def, env): undefined subclass
## "ndiMatrix" of class "replValueSp"; definition not updated
```

```r
library(microbiome)
```

```
## 
## microbiome R package (microbiome.github.com)
##     
## 
## 
##  Copyright (C) 2011-2022 Leo Lahti, 
##     Sudarshan Shetty et al. <microbiome.github.io>
```

```
## 
## Attaching package: 'microbiome'
```

```
## The following object is masked from 'package:ggplot2':
## 
##     alpha
```

```
## The following object is masked from 'package:base':
## 
##     transform
```

```r
library(ape)
```

```
## 
## Attaching package: 'ape'
```

```
## The following object is masked from 'package:dplyr':
## 
##     where
```

```r
metadata <- read_csv("data/qiime2/qiime2_tidy_metadata.csv",show_col_types = FALSE) %>% 
  tibble::column_to_rownames("sample_id") %>% 
  sample_data(metadata)

otutable <- read_csv("data/qiime2/qiime2_tidy_otutable.csv",show_col_types = FALSE) %>% 
  pivot_wider(id_cols = sample_id, names_from = feature, values_from = count) %>% 
  tibble::column_to_rownames("sample_id") %>% 
  otu_table(otutable, taxa_are_rows = FALSE)

taxonomy <- read_csv("data/qiime2/qiime2_tidy_taxonomy.csv", show_col_types = FALSE) %>%
  tibble::column_to_rownames("feature") %>%
  as.matrix() %>% 
  tax_table(taxonomy)

ps_raw_basic <- merge_phyloseq(metadata, otutable, taxonomy)

library(ape)
ps_tree = rtree(ntaxa(ps_raw_basic), rooted=TRUE, tip.label=taxa_names(ps_raw_basic))

ps_raw <- phyloseq::merge_phyloseq(ps_raw_basic, ps_tree)

ps_rel <- phyloseq::transform_sample_counts(ps_raw, function(x){x / sum(x)})

ps_df_raw <- psmelt(ps_raw)

ps_df_rel <- psmelt(ps_raw)

save(ps_tree, ps_raw, ps_rel, ps_df_raw, ps_df_rel, file = "data/qiime2/qiime2_phyloseq_objects.rda")
```

## Review the phyloseq object


```r
ps_raw
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 95 taxa and 128 samples ]
## sample_data() Sample Data:       [ 128 samples by 12 sample variables ]
## tax_table()   Taxonomy Table:    [ 95 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 95 tips and 94 internal nodes ]
```

...
