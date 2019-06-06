#!/usr/bin/env python

rule presentation_six_spt6_western:
    input:
        blot_path = FIGURES["presentation"]["spt6_western"]["blot_path"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt"
    output:
        pdf = "figures/presentation/presentation_six_spt6_western.pdf",
    params:
        height = eval(str(FIGURES["presentation"]["spt6_western"]["height"])),
        width = eval(str(FIGURES["presentation"]["spt6_western"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_six_spt6_western.R"

rule presentation_six_gene_diagram:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_spec"]
    output:
        pdf = "figures/presentation/presentation_six_gene_diagram.pdf",
    params:
        width = eval(str(FIGURES["presentation"]["gene_diagram"]["width"])),
        height = eval(str(FIGURES["presentation"]["gene_diagram"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_six_gene_diagram.R"

