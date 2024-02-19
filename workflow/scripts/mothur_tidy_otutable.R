library(tidyverse, suppressPackageStartupMessages())

# MOTHUR  otutable
mothur_otutable_wider <- read_tsv("data/mothur/mothur_otutable.tsv", show_col_types = FALSE) %>%
  select(-label, -numOtus) %>%
  rename(sample_id = Group)


mothur_tidy_otutable <- mothur_otutable_wider %>% 
  pivot_longer(-sample_id, names_to="otu", values_to="count") %>% 
  filter(count != 0)

write_csv(mothur_tidy_otutable, "data/mothur/mothur_tidy_otutable.csv")
