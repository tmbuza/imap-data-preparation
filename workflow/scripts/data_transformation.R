# Data Tranformation

## Import libraries and saved objects
set.seed(110912)
source("workflow/scripts/common.R")
library(tidyverse, suppressPackageStartupMessages())
cat("\nPhyloseq objects\n\n")

load("data/qiime2_phyloseq_objects.rda", verbose = T)

library(phyloseq)
library(microbiome)
ps_raw
ps_rel


##  No Transformation, Similar to raw abundance
ps_identity <- microbiome::transform(ps_raw, 'identity')


## Relative abundance
ps_rel = phyloseq::transform_sample_counts(ps_raw, function(x){x / sum(x)})


## Arc sine (asin) transformation
x = otu_table(ps_rel)
y <- x/max(x)
ps_asin <- round(asin(sqrt(y)), 6)
ps_asin <- as.matrix(ps_asin)


## Compositional Version, Similar to relative abundance
ps_compositional <- microbiome::transform(ps_raw, 'compositional')


# ## Z-transform for OTUs
# ps_z_otu <- microbiome::transform(ps_raw, 'Z', 'OTU')
# 
# 
# ## Z-transform for Samples
# ps_z_sample <- microbiome::transform(ps_raw, 'Z', 'sample')


# ## Log10 Transform
# ps_log10 <- microbiome::transform(ps_raw, 'log10')


# ## Log10p Transform
# ps_log10p <- microbiome::transform(ps_raw, 'log10p')


# # CLR Transform: Note that small pseudocount is added if data contains zeroes
# ps_clr <- microbiome::transform(ps_raw, 'clr')


## Shift the baseline
ps_shift <- microbiome::transform(ps_raw, 'shift', shift=1)


## Data Scaling
ps_scale <- microbiome::transform(ps_raw, 'scale', scale=1)

## Transformed Objects
save(
  ps_asin, 
  ps_identity,
  ps_compositional, 
  # ps_z_otu, 
  # ps_z_sample, 
  # ps_log10, 
  # ps_log10p, 
  # ps_clr, 
  ps_shift, 
  ps_scale, 
  file = "data/transformed_data.rda")

# load("data/transformed_data.rda", verbose = T)