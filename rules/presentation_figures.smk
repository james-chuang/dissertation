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

rule presentation_six_tfiib_spreading_swc4:
    input:
        tfiib_data = FIGURES["six"]["tfiib_spreading_swc4"]["tfiib_data"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        frame_1 = "figures/presentation/presentation_six_tfiib_spreading_swc0001.svg",
        frame_2 = "figures/presentation/presentation_six_tfiib_spreading_swc0002.svg",
    params:
        height = eval(str(FIGURES["six"]["tfiib_spreading_swc4"]["beamer_height"])),
        width = eval(str(FIGURES["six"]["tfiib_spreading_swc4"]["beamer_width"])),
    conda:
        "../envs/gganimate.yaml"
    script:
        "../scripts/presentation_six_tfiib_spreading_swc4.R"

rule presentation_six_tss_vs_tfiib:
    input:
        genic = FIGURES["six"]["tss_vs_tfiib"]["genic"],
        intragenic = FIGURES["six"]["tss_vs_tfiib"]["intragenic"],
        antisense = FIGURES["six"]["tss_vs_tfiib"]["antisense"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_six_tss_v_tfiib.pdf",
    params:
        height = eval(str(FIGURES["six"]["tss_vs_tfiib"]["beamer_height"])),
        width = eval(str(FIGURES["six"]["tss_vs_tfiib"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_six_tss_vs_tfiib.R"

rule presentation_six_mnase_heatmaps:
    input:
        netseq_data = FIGURES["six"]["mnase_heatmaps"]["netseq_data"],
        mnase_data = FIGURES["six"]["mnase_heatmaps"]["mnase_data"],
        # quant_data = FIGURES["six"]["mnase_heatmaps"]["quant_data"],
        annotation = FIGURES["six"]["mnase_heatmaps"]["annotation"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_six_mnase_heatmaps.pdf",
    params:
        height = eval(str(FIGURES["six"]["mnase_heatmaps"]["beamer_height"])),
        width = eval(str(FIGURES["six"]["mnase_heatmaps"]["beamer_width"])),
        assay = "NET-seq"
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_six_mnase_heatmaps.R"

rule presentation_six_mnase_metagene:
    input:
        mnase_data = FIGURES["six"]["mnase_metagene"]["mnase_data"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        frame_1 = "figures/presentation/presentation_six_mnase_metagene0001.svg",
        frame_2 = "figures/presentation/presentation_six_mnase_metagene0002.svg",
    params:
        height = eval(str(FIGURES["six"]["mnase_metagene"]["beamer_height"])),
        width = eval(str(FIGURES["six"]["mnase_metagene"]["beamer_width"])),
    conda:
        "../envs/gganimate.yaml"
    script:
        "../scripts/presentation_six_mnase_metagene.R"

rule presentation_six_intragenic_mnase_metagenes:
    input:
        mnase_data = FIGURES["six"]["intragenic_mnase_metagenes"]["mnase_data"],
        # gc_data = FIGURES["six"]["intragenic_mnase_metagenes"]["gc_data"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_six_intragenic_mnase_metagenes.pdf",
    params:
        height = eval(str(FIGURES["six"]["intragenic_mnase_metagenes"]["beamer_height"])),
        width = eval(str(FIGURES["six"]["intragenic_mnase_metagenes"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_six_intragenic_mnase_metagenes.R"

rule presentation_six_tss_seqlogos:
    input:
        data_paths = FIGURES["six"]["tss_seqlogos"]["data_paths"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_six_tss_seqlogos.pdf",
    params:
        height = eval(str(FIGURES["six"]["tss_seqlogos"]["beamer_height"])),
        width = eval(str(FIGURES["six"]["tss_seqlogos"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_six_tss_seqlogos.R"

rule presentation_six_intragenic_tata:
    input:
        tata_genic_path = FIGURES["six"]["intragenic_tata"]["tata_genic"],
        tata_intra_path = FIGURES["six"]["intragenic_tata"]["tata_intragenic"],
        tata_anti_path = FIGURES["six"]["intragenic_tata"]["tata_antisense"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        frame_1 = "figures/presentation/presentation_six_intragenic_tata0001.svg",
        frame_2 = "figures/presentation/presentation_six_intragenic_tata0002.svg",
    params:
        height = eval(str(FIGURES["six"]["intragenic_tata"]["beamer_height"])),
        width = eval(str(FIGURES["six"]["intragenic_tata"]["beamer_width"])),
    conda:
        "../envs/gganimate.yaml"
    script:
        "../scripts/presentation_six_intragenic_tata.R"

rule presentation_five_spt5_depletion:
    input:
        spt5_data = FIGURES["five"]["summary_metagenes"]["spt5_data"],
        annotation = FIGURES["five"]["summary_metagenes"]["annotation"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_five_spt5_depletion.pdf",
    params:
        height = eval(str(FIGURES["five"]["spt5_depletion"]["beamer_height"])),
        width = eval(str(FIGURES["five"]["spt5_depletion"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_five_spt5_depletion.R"

rule presentation_five_netseq_meta:
    input:
        netseq_data = FIGURES["five"]["summary_metagenes"]["netseq_data"],
        annotation = FIGURES["five"]["summary_metagenes"]["annotation"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_five_netseq_meta.pdf",
    params:
        height = eval(str(FIGURES["five"]["netseq_meta"]["beamer_height"])),
        width = eval(str(FIGURES["five"]["netseq_meta"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_five_netseq_meta.R"

rule presentation_five_rnaseq_metagene:
    input:
        data_path = FIGURES["five"]["rnaseq_metagene"]["data_path"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_five_rnaseq_metagene.pdf",
    params:
        height = eval(str(FIGURES["five"]["rnaseq_metagene"]["beamer_height"])),
        width = eval(str(FIGURES["five"]["rnaseq_metagene"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_five_rnaseq_metagene.R"

rule convert_svg_to_pdf:
    input:
        svg = "figures/presentation/presentation_{figure}{frame}.svg",
    output:
        pdf = "figures/presentation/presentation_{figure}{frame}.pdf",
    params:
        h_0 = eval(str(FIGURES["six"]["tfiib_spreading_swc4"]["beamer_height"])),
        w_0 = eval(str(FIGURES["six"]["tfiib_spreading_swc4"]["beamer_width"])),
        h_1 = eval(str(FIGURES["six"]["mnase_metagene"]["beamer_height"])),
        w_1 = eval(str(FIGURES["six"]["mnase_metagene"]["beamer_width"])),
    conda:
        "../envs/librsvg.yaml"
    shell: """
        rsvg-convert -f pdf -a -o {output.pdf} {input.svg}
        """

