# QIIME2 composite
library(tidyverse, suppressPackageStartupMessages())

qiime2_tidy_metadata <- read_csv("data/qiime2/qiime2_tidy_metadata.csv", show_col_types = FALSE)
qiime2_tidy_otutable <- read_csv("data/qiime2/qiime2_tidy_features.csv", show_col_types = FALSE)
qiime2_tidy_taxonomy <- read_csv("data/qiime2/qiime2_tidy_taxonomy.csv", show_col_types = FALSE)

qiime2_composite <- inner_join(qiime2_tidy_metadata, qiime2_tidy_otutable, by = "sample_id") %>% 
  inner_join(., qiime2_tidy_taxonomy, by = "feature") %>% 
  group_by(sample_id) %>% 
  mutate(rel_abund = count/sum(count)) %>% 
  ungroup() %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  relocate(count, .before = rel_abund) 

write_csv(qiime2_composite, "data/qiime2/qiime2_composite.csv")