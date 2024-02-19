rule qiime2_phyloseq_objects:
    input:
        metadata=rules.process_qiime2_tidy_metadata.output,
        features=rules.process_qiime2_tidy_features.output,
        taxonomy=rules.process_qiime2_tidy_taxonomy.output,
    output:
        "data/qiime2/qiime2_phyloseq_objects.rda",
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/qiime2_phyloseq_objects.R"




        