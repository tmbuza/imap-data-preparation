# (PART) QIIME2 OUTPUT {-}

# Prepare Qiime2 Output




## Data directories

```r
if (!dir.exists('data')) {dir.create('data')}
if (!dir.exists('data/qiime2')) {dir.create('data/qiime2')}

library(tidyverse, suppressPackageStartupMessages())
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

