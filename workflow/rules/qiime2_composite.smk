# qiime2_metadata
rule qiime2_tidy_metadata:
    input:
        metadata=rules.import_qiime2_metadata.output
    output:
        "data/qiime2/qiime2_tidy_metadata.csv"
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/qiime2_tidy_metadata.R"


# qiime2_features
rule qiime2_tidy_features:
    input:
        features=rules.import_qiime2_features.output
    output:
        "data/qiime2/qiime2_tidy_features.csv"
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/qiime2_tidy_features.R"


# qiime2_taxonomy
rule qiime2_tidy_taxonomy:
    input:
        taxonomy=rules.import_qiime2_taxonomy.output
    output:
        "data/qiime2/qiime2_tidy_taxonomy.csv"
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/qiime2_tidy_taxonomy.R"


# qiime2_composite
rule qiime2_composite_Robject:
    input:
        metadata=rules.qiime2_tidy_metadata.output,
        features=rules.qiime2_tidy_features.output,
        taxonomy=rules.qiime2_tidy_taxonomy.output,
    output:
        "data/qiime2/qiime2_composite.csv"
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/qiime2_composite.R"