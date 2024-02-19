library(tidyverse)
library(vegan)
set.seed(2023)
shared <- read_tsv("~/Dropbox/CDILLC/GIT_REPOS/smda-end2end/data/final.tx.shared", 
                   show_col_types = F) %>%
  select(sample_id = Group, starts_with("Phylo")) %>% 
  pivot_longer(-sample_id, names_to = "otu", values_to = "count") %>%
  group_by(sample_id) %>%
  mutate(total = sum(count)) %>%
  filter(total > 0) %>%
  group_by(otu) %>%
  mutate(total = sum(count)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total)

## Random sampling of pooled data
rand <- shared %>%
  uncount(count) %>%
  mutate(otu = sample(otu)) %>%
  count(sample_id, otu, name="count")

richness <- function(x){
  sum(x>0)
}
shannon <- function(x){
  rabund <- x[x>0]/sum(x) 
  -sum(rabund * log(rabund))
}
simpson <- function(x){
  n <- sum(x)
  sum(x * (x-1) / (n * (n-1)))
}
rarefy <- function(x, sample){
  x <- x[x>0]
  sum(1-exp(lchoose(sum(x) - x, sample) - lchoose(sum(x), sample)))
  
}