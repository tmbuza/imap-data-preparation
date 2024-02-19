# mothur_metadata
rule import_mothur_metadata:
    output:
        "data/mothur/mothur_sample_metadata.tsv"
    conda:
        "../envs/environment.yml"
    shell:
        "bash workflow/scripts/import_mothur_metadata.sh"


# mothur_shared_file
rule import_mothur_otutable:
    output:
        "data/mothur/mothur_otutable.tsv"
    conda:
        "../envs/environment.yml"
    shell:
        "bash workflow/scripts/import_mothur_otutable.sh"


# mothur_taxonomy
rule import_mothur_taxonomy:
    output:
        "data/mothur/mothur_taxonomy.tsv"
    conda:
        "../envs/environment.yml"
    shell:
        "bash workflow/scripts/import_mothur_taxonomy.sh"