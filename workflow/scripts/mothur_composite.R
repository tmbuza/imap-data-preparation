# MOTHUR composite
library(tidyverse, suppressPackageStartupMessages())

mothur_tidy_metadata <- read_csv("data/mothur/mothur_tidy_metadata.csv", show_col_types = FALSE)
mothur_tidy_otutable <- read_csv("data/mothur/mothur_tidy_otutable.csv", show_col_types = FALSE)
mothur_tidy_taxonomy <- read_csv("data/mothur/mothur_tidy_taxonomy.csv", show_col_types = FALSE)

mothur_composite <- inner_join(mothur_tidy_metadata, mothur_tidy_otutable, by = "sample_id") %>% 
  inner_join(., mothur_tidy_taxonomy, by = "otu") %>% 
  group_by(sample_id) %>% 
  mutate(rel_abund = count/sum(count)) %>% 
  ungroup() %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  relocate(count, .before = rel_abund)

write_csv(mothur_composite, "data/mothur/mothur_composite.csv")