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

rule tss_seq_heatmaps:
    input:
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"],
        heatmap_scripts = "scripts/plot_heatmap.R",
        annotation = FIGURES["six"]["tss_seq_heatmaps"]["annotation"],
        tss_sense = FIGURES["six"]["tss_seq_heatmaps"]["tss_sense"],
        tss_antisense = FIGURES["six"]["tss_seq_heatmaps"]["tss_antisense"]
    output:
        pdf = "figures/six/six_tss_seq_heatmaps.pdf"
    params:
        height = eval(str(FIGURES["six"]["tss_seq_heatmaps"]["height"])),
        width = eval(str(FIGURES["six"]["tss_seq_heatmaps"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_tss_seq_heatmaps.R"

rule tfiib_heatmap:
    input:
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"],
        heatmap_scripts = "scripts/plot_heatmap.R",
        annotation = FIGURES["six"]["tfiib_heatmap"]["annotation"],
        tfiib_data = FIGURES["six"]["tfiib_heatmap"]["tfiib_data"],
    output:
        pdf = "figures/six/six_tfiib_heatmap.pdf"
    params:
        height = eval(str(FIGURES["six"]["tfiib_heatmap"]["height"])),
        width = eval(str(FIGURES["six"]["tfiib_heatmap"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_tfiib_heatmap.R"

rule tss_diffexp_summary:
    input:
        in_genic = FIGURES["six"]["tss_diffexp_summary"]["genic"],
        in_intra = FIGURES["six"]["tss_diffexp_summary"]["intragenic"],
        in_anti = FIGURES["six"]["tss_diffexp_summary"]["antisense"],
        in_inter = FIGURES["six"]["tss_diffexp_summary"]["intergenic"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_tss_diffexp_summary.pdf",
    params:
        fdr_cutoff_tss = config["tss_fdr_cutoff"],
        height = eval(str(FIGURES["six"]["tss_diffexp_summary"]["height"])),
        width = eval(str(FIGURES["six"]["tss_diffexp_summary"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_tss_diffexp_summary.R"

rule tss_expression_levels:
    input:
        tss_genic = FIGURES["six"]["tss_expression_levels"]["tss_genic"],
        tss_intragenic = FIGURES["six"]["tss_expression_levels"]["tss_intragenic"],
        tss_antisense = FIGURES["six"]["tss_expression_levels"]["tss_antisense"],
        tss_intergenic = FIGURES["six"]["tss_expression_levels"]["tss_intergenic"],
        # tfiib_genic = FIGURES["figure_one"]["one_d"]["tfiib_genic"],
        # tfiib_intragenic = FIGURES["figure_one"]["one_d"]["tfiib_intragenic"],
        # tfiib_intergenic = FIGURES["figure_one"]["one_d"]["tfiib_intergenic"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_tss_expression_levels.pdf",
    params:
        height = eval(str(FIGURES["six"]["tss_expression_levels"]["height"])),
        width = eval(str(FIGURES["six"]["tss_expression_levels"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_tss_expression_levels.R"


