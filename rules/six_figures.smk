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

rule aat_assay_comparison:
    input:
        data_path = FIGURES["six"]["aat_assay_comparison"]["data_path"],
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"]
    output:
        pdf = "figures/six/six_aat_assay_comparison.pdf",
    params:
        width = eval(str(FIGURES["six"]["aat_assay_comparison"]["width"])),
        height = eval(str(FIGURES["six"]["aat_assay_comparison"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_aat2_assay_comparisons.R"

