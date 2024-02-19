# qiime2_sample_metadata
rule import_qiime2_metadata:
    output:
        "data/qiime2/qiime2_sample_metadata.tsv"
    conda:
        "../envs/environment.yml"
    shell:
        "bash workflow/scripts/import_qiime2_metadata.sh"


# qiime2_features
rule import_qiime2_features:
    output:
        "data/qiime2/qiime2_features.tsv"
    conda:
        "../envs/environment.yml"
    shell:
        "bash workflow/scripts/import_qiime2_features.sh"


# qiime2_taxonomy
rule import_qiime2_taxonomy:
    output:
        "data/qiime2/qiime2_taxonomy.tsv"
    conda:
        "../envs/environment.yml"
    shell:
        "bash workflow/scripts/import_qiime2_taxonomy.sh"

