# Import biom file to phyloseq

library(biomformat)

load("data/qiime2_phyloseq_objects.rda", verbose = TRUE)
ps_biom <- read_biom("../imap-bioinformatics-qiime2/qiime2_process/transformed/feature-table.biom")
