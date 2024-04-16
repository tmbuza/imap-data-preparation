# (PART) MOTHUR OUTPUT {-}

# Prepare Mothur Output



## Data directories

```r
if (!dir.exists('data')) {dir.create('data')}
if (!dir.exists('data/mothur')) {dir.create('data/mothur')}

library(tidyverse, suppressPackageStartupMessages())
```

## Mothur  metadata

```r
read_tsv("data/mothur/mothur_sample_metadata.tsv", show_col_types = FALSE) %>% 
  write_csv("data/mothur/mothur_tidy_metadata.csv")
```


## Mothur  otutable

```r
read_tsv("../imap-mothur-bioinformatics/mothur_process/asv_analysis/final.asv.shared", skip = 0, show_col_types = FALSE) %>%
  dplyr::rename(group="Group") %>% 
  dplyr::select(-c(label, numASVs)) %>% 
  mutate_at(2:ncol(.), as.numeric) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  pivot_longer(-group, names_to = "OTU", values_to = "count") %>% 
  write_csv("data/mothur/mothur_tidy_otutable.csv")
```


## Mothur taxonomy

```r
read_tsv("../imap-mothur-bioinformatics/mothur_process/asv_analysis/final.asv.ASV.cons.taxonomy", show_col_types=FALSE) %>% 
  distinct() %>%
  dplyr::select(-Size) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  mutate(Taxonomy = gsub("\\(100\\)", "", Taxonomy)) %>%  
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";") %>% 
  write_csv("data/mothur/mothur_tidy_taxonomy.csv")
```


## Mothur composite

```r
# Mothur composite
library(tidyverse, suppressPackageStartupMessages())

mothur_tidy_metadata <- read_csv("data/mothur/mothur_tidy_metadata.csv", show_col_types = FALSE)
mothur_tidy_otutable <- read_csv("data/mothur/mothur_tidy_otutable.csv", show_col_types = FALSE)
mothur_tidy_taxonomy <- read_csv("data/mothur/mothur_tidy_taxonomy.csv", show_col_types = FALSE)

mothur_composite <- inner_join(mothur_tidy_metadata, mothur_tidy_otutable, by = "group") %>% 
  inner_join(., mothur_tidy_taxonomy, by = "OTU") %>% 
  group_by(group) %>% 
  mutate(rel_abund = count/sum(count)) %>% 
  ungroup() %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  relocate(count, .before = rel_abund) 

write_csv(mothur_composite, "data/mothur/mothur_composite.csv")
```

