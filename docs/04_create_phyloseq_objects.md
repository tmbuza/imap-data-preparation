# (PART) PHYLOSEQ OBJECTS {-}

# Creating phyloseq Objects

In this section, we demonstrate how to create phyloseq objects to organize and structure microbiome data for analysis. These objects serve as the foundation for downstream analyses and visualization.



## From mothur output

```r
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ape)

metadata <- read_csv("data/mothur/mothur_tidy_metadata.csv",show_col_types = FALSE) %>% 
  tibble::column_to_rownames("group") %>% 
  sample_data(metadata)

otutable <- read_csv("data/mothur/mothur_tidy_otutable.csv",show_col_types = FALSE) %>% 
  pivot_wider(id_cols = group, names_from = OTU, values_from = count) %>% 
  tibble::column_to_rownames("group") %>% 
  otu_table(otutable, taxa_are_rows = FALSE)

taxonomy <- read_csv("data/mothur/mothur_tidy_taxonomy.csv", show_col_types = FALSE) %>%
  tibble::column_to_rownames("OTU") %>%
  as.matrix() %>% 
  tax_table(taxonomy)

mt_ps_raw_basic <- merge_phyloseq(metadata, otutable, taxonomy)

library(ape)
mt_ps_tree = rtree(ntaxa(mt_ps_raw_basic), rooted=TRUE, tip.label=taxa_names(mt_ps_raw_basic))

mt_ps_raw <- phyloseq::merge_phyloseq(mt_ps_raw_basic, mt_ps_tree)

mt_ps_rel <- phyloseq::transform_sample_counts(mt_ps_raw, function(x){x / sum(x)})

mt_ps_df_raw <- psmelt(mt_ps_raw)

mt_ps_df_rel <- psmelt(mt_ps_raw)

save(mt_ps_tree, mt_ps_raw, mt_ps_rel, mt_ps_df_raw, mt_ps_df_rel, file = "data/mothur_phyloseq_objects.rda")
```


### Review the Mothur phyloseq object


```r
mt_ps_raw
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 920 taxa and 10 samples ]
sample_data() Sample Data:       [ 10 samples by 2 sample variables ]
tax_table()   Taxonomy Table:    [ 920 taxa by 6 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 920 tips and 919 internal nodes ]
```

## From QIIME2 output

```r
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ape)

metadata <- read_csv("data/qiime2/qiime2_tidy_metadata.csv",show_col_types = FALSE) %>% 
  tibble::column_to_rownames("sample_id") %>% 
  sample_data(metadata)

otutable <- read_csv("data/qiime2/qiime2_tidy_otutable.csv",show_col_types = FALSE) %>% 
  pivot_wider(id_cols = sample_id, names_from = feature, values_from = count) %>% 
  tibble::column_to_rownames("sample_id") %>% 
  otu_table(otutable, taxa_are_rows = FALSE)

taxonomy <- read_csv("data/qiime2/qiime2_tidy_taxonomy.csv", show_col_types = FALSE) %>%
  tibble::column_to_rownames("feature") %>%
  as.matrix() %>% 
  tax_table(taxonomy)

q2_ps_raw_basic <- merge_phyloseq(metadata, otutable, taxonomy)

library(ape)
q2_ps_tree = rtree(ntaxa(q2_ps_raw_basic), rooted=TRUE, tip.label=taxa_names(q2_ps_raw_basic))

q2_ps_raw <- phyloseq::merge_phyloseq(q2_ps_raw_basic, q2_ps_tree)

q2_ps_rel <- phyloseq::transform_sample_counts(q2_ps_raw, function(x){x / sum(x)})

q2_ps_df_raw <- psmelt(q2_ps_raw)

q2_ps_df_rel <- psmelt(q2_ps_raw)

save(q2_ps_tree, q2_ps_raw, q2_ps_rel, q2_ps_df_raw, q2_ps_df_rel, file = "data/qiime2_phyloseq_objects.rda")
```


## Review the QIIME2 phyloseq object


```r
q2_ps_raw
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 95 taxa and 128 samples ]
sample_data() Sample Data:       [ 128 samples by 12 sample variables ]
tax_table()   Taxonomy Table:    [ 95 taxa by 6 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 95 tips and 94 internal nodes ]
```



# (PART) EXTERNAL DATASETS {-}

# Phyloseq Objects Demo Datasets
In this section, we utilize publicly available real-world datasets to demonstrate the functionality and analysis workflows of phyloseq objects. These datasets serve as valuable resources for researchers seeking to reproduce and validate their analyses using similar microbiome data. By leveraging these demo datasets, we aim to enhance reproducibility and accessibility in microbiome research and analysis.




## The GlobalPatterns dataset
The GlobalPatterns dataset, sourced from the Earth Microbiome Project (EMP), serves as a comprehensive repository for studying microbial communities worldwide. Here's a breakdown of its key attributes:

- **Source:** GlobalPatterns originates from the Earth Microbiome Project (EMP), collecting samples worldwide.

- **Composition:** It comprises high-throughput sequencing data, revealing the taxonomic composition of microbial communities.

- **Scope:** Samples represent diverse global ecosystems, offering a comprehensive view of microbial biodiversity.

- **Format:** Presented as a phyloseq object in R, it integrates sample metadata and taxonomic abundance for analysis.

- **Utility:** Researchers utilize it for community profiling, abundance testing, and ecological modeling, enhancing understanding of global microbial diversity and function.


```r
library(phyloseq) # for GlobalPatterns dataset
data("GlobalPatterns")

ps_GlobalPatterns <-GlobalPatterns
df_GlobalPatterns <-GlobalPatterns %>% 
  phyloseq::psmelt() %>% 
  tibble::rownames_to_column("sample_id") %>% 
  rename_all(tolower)
```

**GlobalPatterns phyloseq-class**


```r
ps_GlobalPatterns
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]
```

**Columns in GlobalPatterns  dataset**


```r
colnames(df_GlobalPatterns)
 [1] "sample_id"                "otu"                     
 [3] "sample"                   "abundance"               
 [5] "x.sampleid"               "primer"                  
 [7] "final_barcode"            "barcode_truncated_plus_t"
 [9] "barcode_full_length"      "sampletype"              
[11] "description"              "kingdom"                 
[13] "phylum"                   "class"                   
[15] "order"                    "family"                  
[17] "genus"                    "species"                 
cat("\n")
```

## The dietswap dataset
The dietswap dataset, available through the microbiome R package, offers valuable insights into the effects of dietary interventions on microbial communities. Here's an overview of its attributes:

- **Source:** The dietswap dataset is derived from research examining the impact of dietary changes on the human gut microbiome.

- **Composition:** It comprises high-throughput sequencing data, providing insights into the taxonomic composition and dynamics of microbial communities in response to dietary alterations.

- **Scope:** Samples are obtained from human participants undergoing dietary interventions, enabling researchers to explore how different diets influence microbial diversity and function within the gut microbiome.

- **Format:** Presented in a format suitable for microbiome analysis, the dataset includes sample metadata and taxonomic abundance data, facilitating comprehensive analyses of microbial community dynamics.

- **Utility:** Researchers utilize the dietswap dataset to investigate the effects of dietary interventions on gut microbiome composition, contributing to a better understanding of the intricate interactions between diet, host physiology, and microbial ecology.


```r
library(microbiome) # for dietswap dataset
data("dietswap")

dietswap
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
```

> Note: dietswap is missing the phylo_tree slot, we can construct it and add it like so:


```r
library(microbiome)
data('dietswap')
ps_raw_basic <- dietswap

library(ape)
ps_tree = rtree(ntaxa(ps_raw_basic), rooted=TRUE, tip.label=taxa_names(ps_raw_basic))
ps_dietswap <- phyloseq::merge_phyloseq(ps_raw_basic, ps_tree)

df_dietswap <-ps_dietswap %>% 
  phyloseq::psmelt() %>% 
  tibble::rownames_to_column("sample_id") %>% 
  dplyr::select(-9) %>% 
  rename_all(tolower)
```


**Dietswap phyloseq-class**


```r
ps_dietswap
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
```

**Columns in Dietswap  dataset**


```r
colnames(df_dietswap)
 [1] "sample_id"              "otu"                    "sample"                
 [4] "abundance"              "subject"                "sex"                   
 [7] "nationality"            "group"                  "timepoint"             
[10] "timepoint.within.group" "bmi_group"              "phylum"                
[13] "family"                 "genus"                 
```

## Caporaso dataset

The Caporaso dataset provides valuable insights into microbial marker genes and their associations with various biological factors.

- **Source:** The Caporaso dataset is included in the `microbiomeMarker` R package and is derived from research led by Dr. J. Gregory Caporaso. It aims to provide insights into microbial marker genes and their associations with various biological factors.

- **Composition:** The dataset comprises high-throughput sequencing data, focusing on microbial marker genes from diverse biological samples. It offers valuable information regarding the taxonomic composition and functional potential of microbial communities.

- **Scope:** Samples are collected from a range of environments, including but not limited to human microbiomes, environmental samples, and animal microbiomes. This diversity allows researchers to explore microbial diversity across different ecosystems and conditions.

- **Format:** The dataset is structured to facilitate microbiome marker analysis, with sample metadata and taxonomic abundance data included. This format enables researchers to conduct comprehensive analyses of microbial marker genes and their associations with environmental or biological factors.

- **Utility:** Researchers utilize the Caporaso dataset to investigate microbial marker genes' roles in various ecosystems, such as host-microbiome interactions, environmental responses, and disease states. The dataset contributes to a better understanding of microbial ecology and its implications for human health and environmental management.


```r
library(microbiomeMarker) # for caporaso dataset
data("caporaso")

ps_caporaso <-caporaso
df_caporaso <-caporaso %>% 
  phyloseq::psmelt() %>% 
  tibble::rownames_to_column("sample_id") %>% 
  rename_all(tolower)
```


**Caporaso phyloseq-class**


```r
ps_caporaso
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 3426 taxa and 34 samples ]
sample_data() Sample Data:       [ 34 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 3426 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 3426 tips and 3424 internal nodes ]
```

**Columns in Caporaso  dataset**


```r
colnames(df_caporaso)
 [1] "sample_id"                "otu"                     
 [3] "sample"                   "abundance"               
 [5] "sampletype"               "year"                    
 [7] "month"                    "day"                     
 [9] "subject"                  "reportedantibioticusage" 
[11] "dayssinceexperimentstart" "description"             
[13] "kingdom"                  "phylum"                  
[15] "class"                    "order"                   
[17] "family"                   "genus"                   
[19] "species"                 
```


## Kostic_CRC dataset

The Kostic_CRC dataset provides valuable insights into the gut microbiome of individuals diagnosed with colorectal cancer (CRC).

- **Source:** The Kostic_CRC dataset is included in the `microbiomeMarker` R package. It is derived from research conducted by the Kostic Lab and focuses on investigating the gut microbiome in colorectal cancer (CRC) patients.

- **Composition:** This dataset comprises high-throughput sequencing data, specifically targeting the gut microbiome of individuals with colorectal cancer. It provides insights into the taxonomic composition and potential functional characteristics of microbial communities associated with CRC.

- **Scope:** Samples are collected from individuals diagnosed with colorectal cancer, allowing researchers to explore the microbial diversity and potential biomarkers associated with CRC development and progression.

- **Format:** The dataset is structured to facilitate microbiome marker analysis, including sample metadata and taxonomic abundance data. This format enables researchers to conduct comprehensive analyses of microbial community dynamics in colorectal cancer.

- **Utility:** Researchers utilize the Kostic_CRC dataset to investigate the role of the gut microbiome in colorectal cancer pathogenesis, prognosis, and treatment response. The dataset contributes to a deeper understanding of the complex interplay between the gut microbiome and colorectal cancer biology, potentially leading to novel diagnostic or therapeutic strategies.



```r
library(microbiomeMarker) # for kostic_crc dataset
data("kostic_crc")

ps_raw_basic <- kostic_crc

library(ape)
ps_tree = rtree(ntaxa(ps_raw_basic), rooted=TRUE, tip.label=taxa_names(ps_raw_basic))
ps_kostic_crc <- phyloseq::merge_phyloseq(ps_raw_basic, ps_tree)

df_kostic_crc <-kostic_crc %>% 
  phyloseq::psmelt() %>% 
  tibble::rownames_to_column("sample_id") %>% 
  rename_all(tolower)
```


**Kostic_crc phyloseq-class**


```r
ps_kostic_crc
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 2505 taxa and 177 samples ]
sample_data() Sample Data:       [ 177 samples by 71 sample variables ]
tax_table()   Taxonomy Table:    [ 2505 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 2505 tips and 2504 internal nodes ]
```

**Columns in Kostic_crc  dataset**


```r
colnames(df_kostic_crc)
 [1] "sample_id"                     "otu"                          
 [3] "sample"                        "abundance"                    
 [5] "x.sampleid"                    "barcodesequence"              
 [7] "linkerprimersequence"          "necrosis_percent"             
 [9] "target_subfragment"            "assigned_from_geo"            
[11] "experiment_center"             "title"                        
[13] "run_prefix"                    "age"                          
[15] "normal_equivalent_percent"     "fibroblast_and_vessel_percent"
[17] "depth"                         "treatment"                    
[19] "age_at_diagnosis"              "common_name"                  
[21] "host_common_name"              "body_site"                    
[23] "elevation"                     "reports_received"             
[25] "cea"                           "pcr_primers"                  
[27] "collection_date"               "altitude"                     
[29] "env_biome"                     "sex"                          
[31] "platform"                      "race"                         
[33] "bsp_diagnosis"                 "study_center"                 
[35] "country"                       "chemotherapy"                 
[37] "year_of_death"                 "ethnicity"                    
[39] "anonymized_name"               "taxon_id"                     
[41] "sample_center"                 "samp_size"                    
[43] "year_of_birth"                 "original_diagnosis"           
[45] "age_unit"                      "study_id"                     
[47] "experiment_design_description" "description_duplicate"        
[49] "diagnosis"                     "body_habitat"                 
[51] "sequencing_meth"               "run_date"                     
[53] "histologic_grade"              "longitude"                    
[55] "env_matter"                    "target_gene"                  
[57] "env_feature"                   "key_seq"                      
[59] "body_product"                  "tumor_percent"                
[61] "library_construction_protocol" "region"                       
[63] "run_center"                    "tumor_type"                   
[65] "bsp_notes"                     "radiation_therapy"            
[67] "inflammation_percent"          "host_subject_id"              
[69] "pc3"                           "latitude"                     
[71] "osh_diagnosis"                 "stage"                        
[73] "primary_disease"               "host_taxid"                   
[75] "description"                   "kingdom"                      
[77] "phylum"                        "class"                        
[79] "order"                         "family"                       
[81] "genus"                         "species"                      
```

## Save objects for transformation and exploration {#save_ps_objects}


```r
save(df_GlobalPatterns, 
     df_dietswap,  
     df_caporaso,
     df_kostic_crc,     
     file = "data/dataframe_objects.rda")

save(ps_GlobalPatterns, 
     ps_dietswap,
     ps_caporaso,
     ps_kostic_crc,
     file = "data/phyloseq_objects.rda")
```

## Confirm saved objects

```r
load("data/phyloseq_objects.rda", verbose = TRUE)
Loading objects:
  ps_GlobalPatterns
  ps_dietswap
  ps_caporaso
  ps_kostic_crc
```


# (PART) PHYLOSEQ FUNCTIONS {-}

# Reviewing Phyloseq Objects
In this section, we demonstrate basic functions for reviewing and exploring the contents of phyloseq objects. By familiarizing ourselves with these functions, we gain insights into the structure, composition, and metadata associated with phyloseq objects, facilitating effective analysis and interpretation of microbiome data.


```r
library(phyloseq)
library(ape)

ps <- ps_dietswap
ps
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
```

## Taxonomy rank names

```r
rank_names(ps)
[1] "Phylum" "Family" "Genus" 
```

## Number of taxa

```r
ntaxa(ps)
[1] 130
```

## Split by taxon name, e.g. Firmicutes

```r
library(metagMisc) # foe phyloseq_sep_tax()

taxa_stats <- phyloseq_sep_tax(ps, TaxRank = "Phylum", drop_NA = FALSE)
cat("\n")
taxa_stats$Firmicutes
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 76 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 76 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 76 tips and 75 internal nodes ]
cat("\n")
taxa_stats$Bacteroidetes
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 16 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 16 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 16 tips and 15 internal nodes ]
cat("\n")
taxa_stats$Actinobacteria
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 8 tips and 7 internal nodes ]
```


## Number of samples

```r
nsamples(ps)
[1] 222
```

## Sample names

```r
head(sample_names(ps), 20)
 [1] "Sample-1"  "Sample-2"  "Sample-3"  "Sample-4"  "Sample-5"  "Sample-6" 
 [7] "Sample-7"  "Sample-8"  "Sample-9"  "Sample-10" "Sample-11" "Sample-12"
[13] "Sample-13" "Sample-14" "Sample-15" "Sample-16" "Sample-17" "Sample-18"
[19] "Sample-19" "Sample-20"
```


## Sample variables

```r
sample_variables(ps)
[1] "subject"                "sex"                    "nationality"           
[4] "group"                  "sample"                 "timepoint"             
[7] "timepoint.within.group" "bmi_group"             
```

## Abundance table

```r
otu_table(ps)[1:5, 1:5]
OTU Table:          [5 taxa and 5 samples]
                     taxa are rows
                                  Sample-1 Sample-2 Sample-3 Sample-4 Sample-5
Peptococcus niger et rel.                2        4       11       10        2
Lachnospira pectinoschiza et rel.       54      141      118       55       67
Leminorella                              3        1        0        0        0
Streptococcus bovis et rel.            112       41       95       27      110
Sutterella wadsworthia et rel.          77       11       44        7       55
```

## Taxonomy table

```r
tax_table(ps)[1:5, ]
Taxonomy Table:     [5 taxa by 3 taxonomic ranks]:
                                  Phylum           Family                    
Peptococcus niger et rel.         "Firmicutes"     "Clostridium cluster IX"  
Lachnospira pectinoschiza et rel. "Firmicutes"     "Clostridium cluster XIVa"
Leminorella                       "Proteobacteria" "Proteobacteria"          
Streptococcus bovis et rel.       "Firmicutes"     "Bacilli"                 
Sutterella wadsworthia et rel.    "Proteobacteria" "Proteobacteria"          
                                  Genus                              
Peptococcus niger et rel.         "Peptococcus niger et rel."        
Lachnospira pectinoschiza et rel. "Lachnospira pectinoschiza et rel."
Leminorella                       "Leminorella"                      
Streptococcus bovis et rel.       "Streptococcus bovis et rel."      
Sutterella wadsworthia et rel.    "Sutterella wadsworthia et rel."   
```

## Phylogenetic tree if available


## Taxa names

```r
taxa_names(ps)[1:10]
 [1] "Peptococcus niger et rel."         "Lachnospira pectinoschiza et rel."
 [3] "Leminorella"                       "Streptococcus bovis et rel."      
 [5] "Sutterella wadsworthia et rel."    "Oxalobacter formigenes et rel."   
 [7] "Megasphaera elsdenii et rel."      "Escherichia coli et rel."         
 [9] "Bacteroides fragilis et rel."      "Peptostreptococcus micros et rel."
```

## Summarize Phyloseq Object

```r

microbiome::summarize_phyloseq(ps)
[[1]]
[1] "1] Min. number of reads = 1776"

[[2]]
[1] "2] Max. number of reads = 28883"

[[3]]
[1] "3] Total number of reads = 2949085"

[[4]]
[1] "4] Average number of reads = 13284.1666666667"

[[5]]
[1] "5] Median number of reads = 13255"

[[6]]
[1] "7] Sparsity = 0.205647955647956"

[[7]]
[1] "6] Any OTU sum to 1 or less? YES"

[[8]]
[1] "8] Number of singletons = 10"

[[9]]
[1] "9] Percent of OTUs that are singletons \n        (i.e. exactly one read detected across all samples)2.30769230769231"

[[10]]
[1] "10] Number of sample variables are: 8"

[[11]]
[1] "subject"                "sex"                    "nationality"           
[4] "group"                  "sample"                 "timepoint"             
[7] "timepoint.within.group" "bmi_group"             
```

## Sort Samples in ascending and descending

```r
cat("In ascending order\n")
In ascending order
head(sort(phyloseq::sample_sums(ps), decreasing = F))
 Sample-56 Sample-195 Sample-196  Sample-42 Sample-164 Sample-147 
      1776       4477       4606       5006       5541       5624 
cat("\nIn descending order\n")

In descending order
head(sort(phyloseq::sample_sums(ps), decreasing = T))
  Sample-3  Sample-11  Sample-12  Sample-54 Sample-208  Sample-20 
     28883      26895      25238      23529      23131      22413 
```

## Drop Samples Below a Threshold

```r
ps1perc0 <- phyloseq::subset_samples(ps, phyloseq::sample_sums(ps) > 1)
ps1perc0
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
```

## Drop Taxa Below a Threshold

```r
pstaxa01 <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 1, ps)
pstaxa01 
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 120 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 120 taxa by 3 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 120 tips and 119 internal nodes ]
```


# (PART) PHYLOSEQ EXTRA OBJECTS {-}

# Create psExtra Objects

After the creation of phyloseq objects, we proceed to create psExtra objects using functionalities provided by the **microViz** package. The psExtra objects extend the capabilities of phyloseq objects by incorporating additional metadata, experimental conditions, or statistical results, enhancing the depth and breadth of analyses performed using microbiome data.

## Prerequisite

Before creating psExtra objects, it's essential to have pre-existing phyloseq objects containing microbiome data and associated metadata. These phyloseq objects serve as the input data for creating psExtra objects. For example, pre-existing phyloseq objects such as ps_dietswap, ps_GlobalPatterns, ps_caporaso, and ps_kostic_crc prepared in the previous sections should be available.

## Load pre-existing phyloseq objects

```r
load("data/phyloseq_objects.rda", verbose = TRUE)
Loading objects:
  ps_GlobalPatterns
  ps_dietswap
  ps_caporaso
  ps_kostic_crc
```

## Aggregate taxonomic data

```r
# Load required packages
library(phyloseq) # For handling microbiome data
library(microViz) # For advanced microbiome visualization

# Assuming ps_dietswap is already loaded or available

# Aggregate taxonomic data
psextra_clr_dietswap <- tax_agg(
  ps = ps_dietswap,
  rank = "Genus") %>% 
  tax_transform(trans = "clr")

# Display a subset of the transformed data
cat("\nCentered log Ratio (CLR) Transformation\n")

Centered log Ratio (CLR) Transformation
psextra_clr_dietswap %>% 
  otu_get() %>%
  .[1:4, 1:3]
OTU Table:          [3 taxa and 4 samples]
                     taxa are columns
         Peptococcus niger et rel. Lachnospira pectinoschiza et rel.
Sample-1               -0.87049643                          2.357387
Sample-2               -0.87803411                          2.599687
Sample-3                0.02051085                          2.353077
Sample-4                0.17534649                          1.861585
         Leminorella
Sample-1  -0.4880348
Sample-2  -2.0404856
Sample-3  -3.1149834
Sample-4  -3.6226887

# Aggregate taxonomic data
psextra_id_dietswap <- tax_agg(
  ps = ps_dietswap,
  rank = "Genus") %>% 
  tax_transform(trans = "identity")

# Display a subset of the transformed data
cat("\nIdentity Transformation\n")

Identity Transformation
psextra_id_dietswap %>% 
  otu_get() %>%
  .[1:4, 1:3]
OTU Table:          [3 taxa and 4 samples]
                     taxa are columns
         Peptococcus niger et rel. Lachnospira pectinoschiza et rel.
Sample-1                         2                                54
Sample-2                         4                               141
Sample-3                        11                               118
Sample-4                        10                                55
         Leminorella
Sample-1           3
Sample-2           1
Sample-3           0
Sample-4           0

# Aggregate taxonomic data
psextra_log10p_dietswap <- tax_agg(
  ps = ps_dietswap,
  rank = "Genus") %>% 
  tax_transform(trans = "log10p")

# Display a subset of the transformed data
cat("\nLog10p Transformation\n")

Log10p Transformation
psextra_log10p_dietswap %>% 
  otu_get() %>%
  .[1:4, 1:3]
OTU Table:          [3 taxa and 4 samples]
                     taxa are columns
         Peptococcus niger et rel. Lachnospira pectinoschiza et rel.
Sample-1                 0.4771213                          1.740363
Sample-2                 0.6989700                          2.152288
Sample-3                 1.0791812                          2.075547
Sample-4                 1.0413927                          1.748188
         Leminorella
Sample-1     0.60206
Sample-2     0.30103
Sample-3     0.00000
Sample-4     0.00000
```


## Confirm the psExtra class

```r
summary(psextra_clr_dietswap)
 Length   Class    Mode 
      1 psExtra      S4 
summary(psextra_id_dietswap)
 Length   Class    Mode 
      1 psExtra      S4 
summary(psextra_log10p_dietswap)
 Length   Class    Mode 
      1 psExtra      S4 
```

## Save psExtra to R object

```r
save(psextra_clr_dietswap, psextra_id_dietswap, psextra_log10p_dietswap, file = "data/phyloseq_extra_objects.rda")
```


# (PART) CLEANING PHYLOSEQ OBJECT {-}

# Pipeline for Cleaning Phyloseq Objects

Before proceeding with downstream analyses, it's a good practice to inspect and clean your phyloseq object to minimize any potential errors or biases. By ensuring data integrity and consistency, you can enhance the reliability of your results.

In this section, we'll outline a pipeline for inspecting and cleaning a phyloseq object. We'll thoroughly examine its contents and rebuild the object as needed to prepare it for subsequent analyses.


## Demo with the `dietswap` dataset


```r
# Set seed for reproducibility
set.seed(110912)

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome) # for dietswap ps object
library(microViz) # for ps_melt function and psExtra object.
library(ape)

# Load the dietswap dataset from the microbiome package
data("dietswap", package = "microbiome")

ps <- ps_dietswap

# Apply string replacement to specified columns
ps_df <- ps_melt(ps)  %>% # melt using microViz function ps_melt()
mutate(Genus = str_replace_all(Genus, " et rel.", ""),
       Genus = str_replace_all(Genus, " at rel.", ""),
       # Species = replace(Species, is.na(Species), 0),
       Genus = replace(Genus, is.na(Genus), 0),
       Family = replace(Family, is.na(Family), 0),
       # Order = replace(Order, is.na(Order), 0),
       # Class = replace(Genus, is.na(Class), 0),
       Phylum = replace(Phylum, is.na(Phylum), 0)
      )

otu_table <- ps_df %>%
pivot_wider(id_cols = "OTU", names_from = "Sample", values_from = "Abundance") %>%
mutate(OTU = paste0("OTU", sprintf("%03d", 1:nrow(otu_table(ps))))) %>% # just for neat OTU names
tibble::column_to_rownames("OTU") %>%
otu_table(otutable, taxa_are_rows = TRUE)

tax_table <- ps_df %>%
select("OTU", "Phylum", "Family", "Genus") %>%
distinct() %>%
mutate(OTU = paste0("OTU", sprintf("%03d", 1:nrow(otu_table(ps))))) %>%
tibble::column_to_rownames("OTU") %>%
as.matrix() %>%
tax_table(taxtable)

# Check if rowname of otutable and that of taxonomy are equal
x1 = otu_table %>%
as.data.frame() %>%
rownames_to_column("OTU")


x2 = tax_table %>%
as.data.frame() %>%
rownames_to_column("OTU")

identical(x1$OTU, x2$OTU)
[1] TRUE


## Creating a new phyloseq object after alteration

ps_basic<- merge_phyloseq(otu_table, tax_table, sample_data(ps))
ps_tree = rtree(ntaxa(ps_basic), rooted=TRUE, tip.label=taxa_names(ps_basic))
ps_raw <- phyloseq::merge_phyloseq(ps_basic, ps_tree)
ps_rel <- phyloseq::transform_sample_counts(ps_raw, function(x){x / sum(x)})

# # Check points
# otu_table(ps_raw)
# tax_table(ps_raw)
# sample_data(ps_raw)
# phy_tree(ps_raw)

## Create a Phylosic extra object

library(microViz)

psextra_raw <- tax_agg(ps = ps_raw, rank = "Genus")

cat("\nClass of ps_raw\n\n")

Class of ps_raw
summary(ps_raw)
  Length    Class     Mode 
       1 phyloseq       S4 

cat("\nClass of psextra_raw\n")

Class of psextra_raw
summary(psextra_raw)
 Length   Class    Mode 
      1 psExtra      S4 

otu_table(psextra_raw)[1:3, 1:3] %>% as.data.frame()
                           Sample-1 Sample-2 Sample-3
Parabacteroides distasonis      187       28       17
Eggerthella lenta                 3       10        8
Coprococcus eutactus             49      355      477

psextra_rel <- phyloseq::transform_sample_counts(psextra_raw, function(x){x / sum(x)})

otu_table(ps_rel)[1:3, 1:3] %>% as.data.frame()
           Sample-1     Sample-2     Sample-3
OTU016 0.0221066320 0.0013278953 0.0005885815
OTU038 0.0003546519 0.0004742483 0.0002769795
OTU055 0.0057926469 0.0168358152 0.0165149050


## Getting clean dataframe from psExtra object
ps_df <- psmelt(psextra_raw) %>%
group_by(Sample) %>%
mutate(total = sum(Abundance)) %>%
filter(total > 0) %>%
filter(Abundance >0) %>%
group_by(OTU) %>%
mutate(total = sum(Abundance)) %>%
filter(total != 0) %>%
ungroup() %>%
select(-total) %>% as.data.frame() %>%
group_by(Sample) %>%
mutate(rel_abund = Abundance/sum(Abundance)) %>%
ungroup() %>%
relocate(Abundance, .before = rel_abund) %>%
rename_all(~ make.unique(tolower(.), sep = "_")) %>%
rename(sample_id = sample,
count = abundance,
bmi = bmi_group,
grp = group,
timewithingrp = timepoint.within.group) %>%
pivot_longer(cols = c("phylum", "family", "genus", "otu"), names_to = "level", values_to = "taxon") %>%
mutate(taxon = str_replace(string = taxon,
pattern = "(.*)",
replacement = "*\\1*"),
taxon = str_replace(string = taxon,
pattern = "\\*(.*)_unclassified\\*",
replacement = "Unclassified<br>*\\1*"),
taxon = str_replace_all(taxon, "_", " "))
```


## Save clean object


```r
save(ps_raw, ps_rel, psextra_raw, psextra_rel, ps_df, file = "data/phyloseq_raw_rel_psextra_df_objects.rda")
```

