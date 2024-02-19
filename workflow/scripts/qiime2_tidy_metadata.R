library(tidyverse, suppressPackageStartupMessages())

# QIIME2  metadata
qiime2_tidy_metadata <- read_tsv("data/qiime2/qiime2_sample_metadata.tsv", show_col_types = FALSE) %>% 
  rename(sample_id="sample-id")

write_csv(qiime2_tidy_metadata, "data/qiime2/qiime2_tidy_metadata.csv")
