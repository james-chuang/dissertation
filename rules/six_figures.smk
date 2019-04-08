#!/usr/bin/env python

rule gene_diagram:
    input:
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"]
    output:
        pdf = "figures/six/six_gene_diagram.pdf",
    params:
        width = eval(str(FIGURES["six"]["gene_diagram"]["width"])),
        height = eval(str(FIGURES["six"]["gene_diagram"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_gene_diagram.R"


