#!/usr/bin/env python

rule five_summary_metagenes:
    input:
        netseq_data = FIGURES["five"]["summary_metagenes"]["netseq_data"],
        rnapii_data = FIGURES["five"]["summary_metagenes"]["rnapii_data"],
        spt5_data = FIGURES["five"]["summary_metagenes"]["spt5_data"],
        annotation = FIGURES["five"]["summary_metagenes"]["annotation"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/five/five_summary_metagenes.pdf",
    params:
        height = eval(str(FIGURES["five"]["summary_metagenes"]["height"])),
        width = eval(str(FIGURES["five"]["summary_metagenes"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_summary_metagenes.R"

rule five_rnapii_phosphomark_enrichment:
    input:
        data_paths = FIGURES["five"]["rnapii_phosphomark_enrichment"]["data_paths"],
        annotation = FIGURES["five"]["rnapii_phosphomark_enrichment"]["annotation"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/five/five_rnapii_phosphomark_enrichment.pdf",
    params:
        height = eval(str(FIGURES["five"]["rnapii_phosphomark_enrichment"]["height"])),
        width = eval(str(FIGURES["five"]["rnapii_phosphomark_enrichment"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_rnapii_phosphomark_enrichment.R"

rule five_rnaseq_heatmaps:
    input:
        rnaseq_data = FIGURES["five"]["rnaseq_heatmaps"]["rnaseq_data"],
        annotation = FIGURES["five"]["rnaseq_heatmaps"]["annotation"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
        heatmap_scripts = "scripts/plot_heatmap.R",
    output:
        pdf = "figures/five/five_rnaseq_heatmaps.pdf",
    params:
        height = eval(str(FIGURES["five"]["rnaseq_heatmaps"]["height"])),
        width = eval(str(FIGURES["five"]["rnaseq_heatmaps"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_rnaseq_heatmaps.R"

