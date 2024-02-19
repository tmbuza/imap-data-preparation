rule plot_taxa_barplot:
    input:
        "data/mothur_composite.csv",
        "data/qiime2_composite.csv",
    output:
        report("figures/taxon_barplot.svg", caption="../report/barplot_1.rst", category="Genus Abundance"),
        "figures/taxon_barplot.png",
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/plot_merged_taxa.R"