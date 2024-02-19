library(tidyverse, suppressPackageStartupMessages())

# MOTHUR taxonomy
mothur_tidy_taxonomy <- read_tsv("data/mothur/mothur_taxonomy.tsv", show_col_types=FALSE) %>% 
  rename_all(tolower) %>%
  select(otu, taxonomy) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";unclassified", "_unclassified"),
         taxonomy = str_replace_all(taxonomy, ";unclassified", ""),
         taxonomy = str_replace_all(taxonomy, ";$", "")) %>% 
  # mutate(pretty_otu = str_replace(string = otu,
  #                                 pattern = "Otu0*",
  #                                 replacement = "OTU ")) %>% 
  # relocate(pretty_otu, .after = "otu") %>% 
  separate(taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";")

write_csv(mothur_tidy_taxonomy, "data/mothur/mothur_tidy_taxonomy.csv")

