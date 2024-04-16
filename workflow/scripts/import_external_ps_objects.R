
library(microbiomeMarker) # For caporaso dataset
library(phyloseq) # Also for dietswap dataset 
library(microViz) # For misc conversion, not in CRAN yet, not installed in the environment
library(metagMisc) # For visualization, not in CRAN yet, not installed in the environment
library(metagenomeSeq)
library(tidyverse)
library(corncob) # For IBD example dataset


# Dataset 1
library(phyloseq)
data("GlobalPatterns")
ps_GlobalPatterns <-GlobalPatterns
df_GlobalPatterns <-GlobalPatterns %>% 
  phyloseq::psmelt() %>% 
  tibble::rownames_to_column("sample_id") %>% 
  rename_all(tolower)

colnames(df_GlobalPatterns) %>% as.data.frame()


# Dataset 2
library(corncob)
data("ibd_phylo")
ps_ibd_phylo <-ibd_phylo
df_ibd_phylo <-ibd_phylo %>% 
  phyloseq::psmelt() %>% 
  tibble::rownames_to_column("sample_id") %>% 
  rename_all(tolower)

colnames(df_ibd_phylo) %>% as.data.frame()


# Dataset 3
library(microbiome)
data("dietswap")
ps_dietswap <-dietswap
df_dietswap <-dietswap %>% 
  phyloseq::psmelt() %>% 
  tibble::rownames_to_column("sample_id") %>% 
  rename_all(tolower)

colnames(df_dietswap) %>% as.data.frame()


# Dataset 4
library(microbiomeMarker)
data("caporaso")
ps_caporaso <-caporaso
df_caporaso <-caporaso %>% 
  phyloseq::psmelt() %>% 
  tibble::rownames_to_column("sample_id") %>% 
  rename_all(tolower)


colnames(df_caporaso) %>% as.data.frame()

# Dataset 5
library(microbiomeMarker)
data("kostic_crc")
ps_kostic_crc <-kostic_crc
df_kostic_crc <-kostic_crc %>% 
  phyloseq::psmelt() %>% 
  tibble::rownames_to_column("sample_id") %>% 
  rename_all(tolower)

colnames(df_kostic_crc) %>% as.data.frame()


if (!dir.exists("data")) {dir.create("data")}
if (!dir.exists("data/external")) {dir.create("data/external")}

#-------------------------------------
#-------------------------------------
save(df_GlobalPatterns, 
     df_ibd_phylo, 
     df_dietswap,  
     df_caporaso,
     df_kostic_crc,

     ps_GlobalPatterns, 
     ps_ibd_phylo, 
     ps_dietswap,
     ps_caporaso,
     ps_kostic_crc,
     file = "data/external/phyloseq_objects.rda")

