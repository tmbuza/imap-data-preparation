# (PART) EXTERNAL DATASETS {-}

# Import Preprocessed Microbiome Data



In this section, we are presenting the process of importing authentic preprocessed data, which will serve as the cornerstone for our subsequent analyses. Although synthetic data certainly has its advantages, it can at times introduce biases or limitations that do not fully capture the intricate diversity inherent in real microbiome samples.

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

colnames(df_GlobalPatterns)
```

```
 [1] "sample_id"                "otu"                     
 [3] "sample"                   "abundance"               
 [5] "x.sampleid"               "primer"                  
 [7] "final_barcode"            "barcode_truncated_plus_t"
 [9] "barcode_full_length"      "sampletype"              
[11] "description"              "kingdom"                 
[13] "phylum"                   "class"                   
[15] "order"                    "family"                  
[17] "genus"                    "species"                 
```

```r
cat("\n")
```

```r
ps_GlobalPatterns
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]
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

ps_dietswap <-dietswap
df_dietswap <-dietswap %>% 
  phyloseq::psmelt() %>% 
  tibble::rownames_to_column("sample_id") %>% 
  dplyr::select(-9) %>% 
  rename_all(tolower)

colnames(df_dietswap)
```

```
 [1] "sample_id"              "otu"                    "sample"                
 [4] "abundance"              "subject"                "sex"                   
 [7] "nationality"            "group"                  "timepoint"             
[10] "timepoint.within.group" "bmi_group"              "phylum"                
[13] "family"                 "genus"                 
```

```r
cat("\n")
```

```r
ps_dietswap
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
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

colnames(df_caporaso)
```

```
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

```r
cat("\n")
```

```r
ps_caporaso
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 3426 taxa and 34 samples ]
sample_data() Sample Data:       [ 34 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 3426 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 3426 tips and 3424 internal nodes ]
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

ps_kostic_crc <-kostic_crc
df_kostic_crc <-kostic_crc %>% 
  phyloseq::psmelt() %>% 
  tibble::rownames_to_column("sample_id") %>% 
  rename_all(tolower)

colnames(df_kostic_crc)
```

```
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

```r
cat("\n")
```

```r
ps_kostic_crc
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 2505 taxa and 177 samples ]
sample_data() Sample Data:       [ 177 samples by 71 sample variables ]
tax_table()   Taxonomy Table:    [ 2505 taxa by 7 taxonomic ranks ]
```


## Save objects for transformation and exploration {#save_ps_objects}


```r
save(df_GlobalPatterns, 
     # df_ibd_phylo, 
     df_dietswap,  
     df_caporaso,
     df_kostic_crc,

     ps_GlobalPatterns, 
     # ps_ibd_phylo, 
     ps_dietswap,
     ps_caporaso,
     ps_kostic_crc,
     file = "data/external/external_ps_objects.rda")
```

## Confirm saved objects

```r
load("data/external/external_ps_objects.rda", verbose = TRUE)
```

```
Loading objects:
  df_GlobalPatterns
  df_dietswap
  df_caporaso
  df_kostic_crc
  ps_GlobalPatterns
  ps_dietswap
  ps_caporaso
  ps_kostic_crc
```

## Review Phyloseq Objects

```r
library(phyloseq)
library(ape)

ps <- ps_dietswap

library(metagMisc)
```

### Taxonomy rank names

```r
rank_names(ps)
```

```
[1] "Phylum" "Family" "Genus" 
```

### Number of taxa

```r
ntaxa(ps)
```

```
[1] 130
```

### Split by taxon name, e.g. Firmicutes

```r
taxa_stats <- phyloseq_sep_tax(ps, TaxRank = "Phylum", drop_NA = FALSE)
cat("\n")
```

```r
taxa_stats$Firmicutes
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 76 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 76 taxa by 3 taxonomic ranks ]
```

```r
cat("\n")
```

```r
taxa_stats$Bacteroidetes
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 16 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 16 taxa by 3 taxonomic ranks ]
```

```r
cat("\n")
```

```r
taxa_stats$Actinobacteria
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8 taxa by 3 taxonomic ranks ]
```


### Number of samples

```r
nsamples(ps)
```

```
[1] 222
```

### Sample names

```r
head(sample_names(ps), 20)
```

```
 [1] "Sample-1"  "Sample-2"  "Sample-3"  "Sample-4"  "Sample-5"  "Sample-6" 
 [7] "Sample-7"  "Sample-8"  "Sample-9"  "Sample-10" "Sample-11" "Sample-12"
[13] "Sample-13" "Sample-14" "Sample-15" "Sample-16" "Sample-17" "Sample-18"
[19] "Sample-19" "Sample-20"
```


### Sample variables

```r
sample_variables(ps)
```

```
[1] "subject"                "sex"                    "nationality"           
[4] "group"                  "sample"                 "timepoint"             
[7] "timepoint.within.group" "bmi_group"             
```

### Abundance table

```r
otu_table(ps)[1:5, 1:5]
```

```
OTU Table:          [5 taxa and 5 samples]
                     taxa are rows
                             Sample-1 Sample-2 Sample-3 Sample-4 Sample-5
Actinomycetaceae                    0        1        0        1        0
Aerococcus                          0        0        0        0        0
Aeromonas                           0        0        0        0        0
Akkermansia                        18       97       67      256       21
Alcaligenes faecalis et rel.        1        2        3        2        2
```

### Taxonomy table

```r
tax_table(ps)[1:5, ]
```

```
Taxonomy Table:     [5 taxa by 3 taxonomic ranks]:
                             Phylum            Family           
Actinomycetaceae             "Actinobacteria"  "Actinobacteria" 
Aerococcus                   "Firmicutes"      "Bacilli"        
Aeromonas                    "Proteobacteria"  "Proteobacteria" 
Akkermansia                  "Verrucomicrobia" "Verrucomicrobia"
Alcaligenes faecalis et rel. "Proteobacteria"  "Proteobacteria" 
                             Genus                         
Actinomycetaceae             "Actinomycetaceae"            
Aerococcus                   "Aerococcus"                  
Aeromonas                    "Aeromonas"                   
Akkermansia                  "Akkermansia"                 
Alcaligenes faecalis et rel. "Alcaligenes faecalis et rel."
```

### Phylogenetic tree if available


### Taxa names

```r
taxa_names(ps)[1:10]
```

```
 [1] "Actinomycetaceae"                  "Aerococcus"                       
 [3] "Aeromonas"                         "Akkermansia"                      
 [5] "Alcaligenes faecalis et rel."      "Allistipes et rel."               
 [7] "Anaerobiospirillum"                "Anaerofustis"                     
 [9] "Anaerostipes caccae et rel."       "Anaerotruncus colihominis et rel."
```

### Summarize Phyloseq Object

```r
microbiome::summarize_phyloseq(ps)
```

```
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


### Sort Samples in ascending and descending

```r
cat("In ascending order\n")
```

```
In ascending order
```

```r
head(sort(phyloseq::sample_sums(ps), decreasing = F))
```

```
 Sample-56 Sample-195 Sample-196  Sample-42 Sample-164 Sample-147 
      1776       4477       4606       5006       5541       5624 
```

```r
cat("\nIn descending order\n")
```

```

In descending order
```

```r
head(sort(phyloseq::sample_sums(ps), decreasing = T))
```

```
  Sample-3  Sample-11  Sample-12  Sample-54 Sample-208  Sample-20 
     28883      26895      25238      23529      23131      22413 
```

### Drop Samples Below a Threshold

```r
ps1perc0 <- phyloseq::subset_samples(ps, phyloseq::sample_sums(ps) > 1)
ps1perc0
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
```

### Drop Taxa Below a Threshold

```r
pstaxa01 <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 1, ps)
pstaxa01 
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 120 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 120 taxa by 3 taxonomic ranks ]
```




