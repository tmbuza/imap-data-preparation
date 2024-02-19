# mothur_metadata
rule mothur_tidy_metadata:
    input:
        metadata=rules.import_mothur_metadata.output
    output:
        "data/mothur/mothur_tidy_metadata.csv"
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/mothur_tidy_metadata.R"


# mothur_shared_file
rule mothur_tidy_otutable:
    input:
        otutable=rules.import_mothur_otutable.output
    output:
        "data/mothur/mothur_tidy_otutable.csv"
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/mothur_tidy_otutable.R"


# mothur_taxonomy
rule mothur_tidy_taxonomy:
    input:
        taxonomy=rules.import_mothur_taxonomy.output
    output:
        "data/mothur/mothur_tidy_taxonomy.csv"
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/mothur_tidy_taxonomy.R"


# mothur_composite
rule mothur_composite_Robject:
    input:
        metadata=rules.mothur_tidy_metadata.output,
        otutable=rules.mothur_tidy_otutable.output,
        taxonomy=rules.mothur_tidy_taxonomy.output,
    output:
        "data/mothur/mothur_composite.csv"
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/mothur_composite.R"