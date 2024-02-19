library(tidyverse, suppressPackageStartupMessages())

# QIIME2 taxonomy
qiime2_tidy_taxonomy <- read_tsv("data/qiime2/qiime2_taxonomy.tsv", show_col_types=FALSE) %>% 
  rename(feature="Feature ID") %>% 
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
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";")

write_csv(qiime2_tidy_taxonomy, "data/qiime2/qiime2_tidy_taxonomy.csv")

