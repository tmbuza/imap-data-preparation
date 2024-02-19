library(tidyverse, suppressPackageStartupMessages())

# MOTHUR  metadata
mothur_tidy_metadata <- read_tsv("data/mothur/mothur_sample_metadata.tsv", show_col_types = FALSE) %>% 
  rename_all(tolower) %>% 
  dplyr::filter(!grepl("Mock*", group))

mothur_tidy_metadata %>% 
  rename(sample_id = 1) %>% 
  select(1, 2) %>% 
write_csv("data/mothur/mothur_tidy_metadata.csv")
