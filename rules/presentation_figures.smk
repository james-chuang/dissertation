#!/usr/bin/env python

rule presentation_six_spt6_western:
    input:
        blot_path = FIGURES["six"]["spt6_western"]["blot_path"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt"
    output:
        pdf = "figures/presentation/presentation_spt6_western.pdf",
    params:
        height = eval(str(FIGURES["six"]["spt6_western"]["beamer_height"])),
        width = eval(str(FIGURES["six"]["spt6_western"]["beamer_width"])),
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
        width = eval(str(FIGURES["six"]["gene_diagram"]["beamer_width"])),
        height = eval(str(FIGURES["six"]["gene_diagram"]["beamer_height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_six_gene_diagram.R"

rule presentation_six_aat_assay_comparison:
    input:
        data_path = FIGURES["six"]["aat_assay_comparison"]["data_path"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_six_aat_assay_comparison.pdf",
    params:
        width = eval(str(FIGURES["six"]["aat_assay_comparison"]["beamer_width"])),
        height = eval(str(FIGURES["six"]["aat_assay_comparison"]["beamer_height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_six_aat2_assay_comparisons.R"

rule presentation_six_tss_seq_heatmaps:
    input:
        theme = config["theme_spec"],
        heatmap_scripts = "scripts/plot_heatmap.R",
        annotation = FIGURES["six"]["tss_seq_heatmaps"]["annotation"],
        tss_sense = FIGURES["six"]["tss_seq_heatmaps"]["tss_sense"],
        tss_antisense = FIGURES["six"]["tss_seq_heatmaps"]["tss_antisense"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_six_tss_seq_heatmaps.pdf"
    params:
        height = eval(str(FIGURES["six"]["tss_seq_heatmaps"]["beamer_height"])),
        width = eval(str(FIGURES["six"]["tss_seq_heatmaps"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_six_tss_seq_heatmaps.R"

rule presentation_six_tss_diffexp_summary:
    input:
        in_genic = FIGURES["six"]["tss_diffexp_summary"]["genic"],
        in_intra = FIGURES["six"]["tss_diffexp_summary"]["intragenic"],
        in_anti = FIGURES["six"]["tss_diffexp_summary"]["antisense"],
        in_inter = FIGURES["six"]["tss_diffexp_summary"]["intergenic"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_six_tss_diffexp_summary.pdf",
    params:
        fdr_cutoff_tss = config["tss_fdr_cutoff"],
        height = eval(str(FIGURES["six"]["tss_diffexp_summary"]["beamer_height"])),
        width = eval(str(FIGURES["six"]["tss_diffexp_summary"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_six_tss_diffexp_summary.R"

rule presentation_six_tss_expression_levels:
    input:
        tss_genic = FIGURES["six"]["tss_expression_levels"]["tss_genic"],
        tss_intragenic = FIGURES["six"]["tss_expression_levels"]["tss_intragenic"],
        tss_antisense = FIGURES["six"]["tss_expression_levels"]["tss_antisense"],
        tss_intergenic = FIGURES["six"]["tss_expression_levels"]["tss_intergenic"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_six_tss_expression_levels.pdf",
    params:
        height = eval(str(FIGURES["six"]["tss_expression_levels"]["beamer_height"])),
        width = eval(str(FIGURES["six"]["tss_expression_levels"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_six_tss_expression_levels.R"

rule presentation_six_tfiib_heatmap:
    input:
        theme = config["theme_spec"],
        heatmap_scripts = "scripts/plot_heatmap.R",
        annotation = FIGURES["six"]["tfiib_heatmap"]["annotation"],
        tfiib_data = FIGURES["six"]["tfiib_heatmap"]["tfiib_data"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_six_tfiib_heatmap.pdf"
    params:
        height = eval(str(FIGURES["six"]["tfiib_heatmap"]["beamer_height"])),
        width = eval(str(FIGURES["six"]["tfiib_heatmap"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_six_tfiib_heatmap.R"

