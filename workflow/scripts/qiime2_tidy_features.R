library(tidyverse, suppressPackageStartupMessages())

# QIIME2  feature table
qiime2_otutable_wider <- read_tsv("data/qiime2/qiime2_features.tsv", skip = 1, show_col_types = FALSE) %>%
  rename(feature='#OTU ID') %>%
  select(-starts_with('Mock')) %>% 
  mutate_at(2:ncol(.), as.numeric) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0))

qiime2_tidy_otutable <- qiime2_otutable_wider %>% 
  pivot_longer(-feature, names_to = "sample_id", values_to = "count") %>% 
  relocate(sample_id, .before = feature)

write_csv(qiime2_tidy_otutable, "data/qiime2/qiime2_tidy_features.csv")
