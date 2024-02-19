rule example_phyloseq_objects:
    output:
        "data/external/example.rda",
    conda:
        "../envs/environment.yml"
    shell:
        "touch {output}"


        