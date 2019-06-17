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

rule presentation_six_intragenic_gc:
    input:
        # mnase_data = FIGURES["six"]["intragenic_mnase_metagenes"]["mnase_data"],
        gc_data = FIGURES["six"]["intragenic_mnase_metagenes"]["gc_data"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_six_intragenic_gc.pdf",
    params:
        height = eval(str(FIGURES["six"]["intragenic_gc"]["beamer_height"])),
        width = eval(str(FIGURES["six"]["intragenic_gc"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_six_intragenic_gc.R"

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

rule presentation_five_rnapii_phosphomark_enrichment:
    input:
        data_paths = FIGURES["five"]["rnapii_phosphomark_enrichment"]["data_paths"],
        annotation = FIGURES["five"]["rnapii_phosphomark_enrichment"]["annotation"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_five_rnapii_phosphomark_enrichment.pdf",
    params:
        height = eval(str(FIGURES["five"]["rnapii_phosphomark_enrichment"]["beamer_height"])),
        width = eval(str(FIGURES["five"]["rnapii_phosphomark_enrichment"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_five_rnapii_phosphomark_enrichment.R"

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

rule presentation_five_rnaseq_heatmaps:
    input:
        rnaseq_data = FIGURES["five"]["rnaseq_heatmaps"]["rnaseq_data"],
        annotation = FIGURES["five"]["rnaseq_heatmaps"]["annotation"],
        theme = config["theme_spec"],
        heatmap_scripts = "scripts/plot_heatmap.R",
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_five_rnaseq_heatmaps.pdf",
    params:
        height = eval(str(FIGURES["five"]["rnaseq_heatmaps"]["beamer_height"])),
        width = eval(str(FIGURES["five"]["rnaseq_heatmaps"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_five_rnaseq_heatmaps.R"

rule presentation_five_antisense_heatmaps:
    input:
        theme = config["theme_spec"],
        tssseq_data = FIGURES["five"]["antisense_heatmaps"]["tssseq_data"],
        rnaseq_data = FIGURES["five"]["antisense_heatmaps"]["rnaseq_data"],
        # netseq_data = FIGURES["five"]["antisense_heatmaps"]["netseq_data"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_five_antisense_heatmaps.pdf"
    params:
        height = eval(str(FIGURES["five"]["antisense_heatmaps"]["beamer_height"])),
        width = eval(str(FIGURES["five"]["antisense_heatmaps"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_five_antisense_heatmaps.R"

rule presentation_five_mnase_metagene:
    input:
        theme = config["theme_spec"],
        mnase_data = FIGURES["five"]["presentation_mnase_metagene"]["mnase_data"],
        tss_peaks = FIGURES["five"]["presentation_mnase_metagene"]["tss_peaks"],
        nuc_dyads = FIGURES["five"]["presentation_mnase_metagene"]["nuc_dyads"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_five_mnase_metagene.pdf"
    params:
        height = eval(str(FIGURES["five"]["presentation_mnase_metagene"]["beamer_height"])),
        width = eval(str(FIGURES["five"]["presentation_mnase_metagene"]["beamer_width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_five_mnase_metagene.R"

rule presentation_stress_tfiib_ridgelines:
    input:
        data_path = FIGURES["stress"]["presentation_tfiib_ridgelines"]["data_path"],
        annotations = FIGURES["stress"]["presentation_tfiib_ridgelines"]["annotations"],
        diamide_genic = FIGURES["stress"]["presentation_tfiib_ridgelines"]["diamide_genic"],
        diamide_intragenic = FIGURES["stress"]["presentation_tfiib_ridgelines"]["diamide_intragenic"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_stress_tfiib_ridgelines.pdf",
    params:
        width = eval(str(FIGURES["stress"]["presentation_tfiib_ridgelines"]["beamer_width"])),
        height = eval(str(FIGURES["stress"]["presentation_tfiib_ridgelines"]["beamer_height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_stress_tfiib_ridgelines.R"

rule presentation_stress_tfiib_coverage:
    input:
        data_paths = [v["paths"] for k,v in  FIGURES["stress"]["tfiib_coverage"]["data"].items()],
        transcript_annotation = FIGURES["stress"]["tfiib_coverage"]["transcript_annotation"],
        orf_annotation = FIGURES["stress"]["tfiib_coverage"]["orf_annotation"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_stress_tfiib_coverage.pdf",
    params:
        gene_ids = FIGURES["stress"]["tfiib_coverage"]["data"].keys(),
        upstream = [v["upstream"] for k,v in FIGURES["stress"]["tfiib_coverage"]["data"].items()],
        downstream = [v["downstream"] for k,v in FIGURES["stress"]["tfiib_coverage"]["data"].items()],
        width = eval(str(FIGURES["stress"]["tfiib_coverage"]["beamer_width"])),
        height = eval(str(FIGURES["stress"]["tfiib_coverage"]["beamer_height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_stress_tfiib_coverage.R"

rule presentation_stress_promoter_tss_expression:
    input:
        genic_path = FIGURES["stress"]["promoter_tss_expression"]["genic_path"],
        intra_path = FIGURES["stress"]["promoter_tss_expression"]["intra_path"],
        counts_path = FIGURES["stress"]["promoter_tss_expression"]["counts_path"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_stress_promoter_tss_expression.pdf",
    params:
        width = eval(str(FIGURES["stress"]["promoter_tss_expression"]["beamer_width"])),
        height = eval(str(FIGURES["stress"]["promoter_tss_expression"]["beamer_height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_stress_promoter_tss_expression.R"

rule presentation_stress_promoter_tss_polyenrichment:
    input:
        matched_peaks_genic = FIGURES["stress"]["promoter_tss_polyenrichment"]["matched_peaks_genic"],
        matched_peaks_intra = FIGURES["stress"]["promoter_tss_polyenrichment"]["matched_peaks_intra"],
        polyenrichment_genic = FIGURES["stress"]["promoter_tss_polyenrichment"]["polyenrichment_genic"],
        polyenrichment_intra = FIGURES["stress"]["promoter_tss_polyenrichment"]["polyenrichment_intra"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_stress_promoter_tss_polyenrichment.pdf",
    params:
        fdr_cutoff_tss = config["tss_fdr_cutoff"],
        width = eval(str(FIGURES["stress"]["promoter_tss_polyenrichment"]["beamer_width"])),
        height = eval(str(FIGURES["stress"]["promoter_tss_polyenrichment"]["beamer_height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_stress_promoter_tss_polyenrichment.R"

rule presentation_stress_dsk2_summary:
    input:
        mnase_chip_paths = FIGURES["stress"]["dsk2_summary"]["mnase_chip_paths"],
        tss_path = FIGURES["stress"]["dsk2_summary"]["tss_path"],
        tfiib_path = FIGURES["stress"]["dsk2_summary"]["tfiib_path"],
        tfiib_diffbind_path = FIGURES["stress"]["dsk2_summary"]["tfiib_diffbind_path"],
        transcript_annotation = FIGURES["stress"]["dsk2_summary"]["transcript_annotation"],
        orf_annotation = FIGURES["stress"]["dsk2_summary"]["orf_annotation"],
        dsk2_pace_annotation = FIGURES["stress"]["dsk2_summary"]["dsk2_pace_annotation"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_stress_dsk2_summary.pdf",
    params:
        gene_ids = FIGURES["stress"]["dsk2_summary"]["gene_ids"],
        up_distances = FIGURES["stress"]["dsk2_summary"]["up_distances"],
        down_distances = FIGURES["stress"]["dsk2_summary"]["down_distances"],
        modification = FIGURES["stress"]["dsk2_summary"]["modification"],
        refpointlabel = FIGURES["stress"]["dsk2_summary"]["refpointlabel"],
        trim_pct = FIGURES["stress"]["dsk2_summary"]["trim_pct"],
        x_label_dist = FIGURES["stress"]["dsk2_summary"]["x_label_dist"],
        width = eval(str(FIGURES["stress"]["dsk2_summary"]["beamer_width"])),
        height = eval(str(FIGURES["stress"]["dsk2_summary"]["beamer_height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_stress_dsk2_summary.R"

rule presentation_stress_dsk2_pace_northern:
    input:
        dsk_blot = FIGURES["stress"]["dsk2_pace_northern"]["dsk_blot"],
        snr_blot = FIGURES["stress"]["dsk2_pace_northern"]["snr_blot"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/presentation/presentation_stress_dsk2_pace_northern.pdf",
    params:
        width = eval(str(FIGURES["stress"]["dsk2_pace_northern"]["beamer_width"])),
        height = eval(str(FIGURES["stress"]["dsk2_pace_northern"]["beamer_height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/presentation_stress_dsk2_pace_northern.R"

rule presentation_stress_diamide_fitnesscomp:
    input:
        data_path = FIGURES["stress"]["diamide_fitnesscomp"]["data_path"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        frame_1 = "figures/presentation/presentation_stress_diamide_fitnesscomp0001.svg",
        frame_2 = "figures/presentation/presentation_stress_diamide_fitnesscomp0002.svg",
        frame_3 = "figures/presentation/presentation_stress_diamide_fitnesscomp0003.svg",
    params:
        width = eval(str(FIGURES["stress"]["diamide_fitnesscomp"]["beamer_width"])),
        height = eval(str(FIGURES["stress"]["diamide_fitnesscomp"]["beamer_height"])),
    conda:
        "../envs/gganimate.yaml"
    script:
        "../scripts/presentation_stress_diamide_fitnesscomp.R"

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

