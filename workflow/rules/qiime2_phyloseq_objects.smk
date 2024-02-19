rule qiime2_phyloseq_objects:
    input:
        metadata=rules.qiime2_tidy_metadata.output,
        features=rules.qiime2_tidy_features.output,
        taxonomy=rules.qiime2_tidy_taxonomy.output,
    output:
        "data/qiime2/qiime2_phyloseq_objects.rda",
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/qiime2_phyloseq_objects.R"




        