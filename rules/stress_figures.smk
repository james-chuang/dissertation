#!/usr/bin/env python

rule gasch_comparison:
    input:
        gasch_one = FIGURES["stress"]["gasch_comparison"]["gasch_one"],
        gasch_three = FIGURES["stress"]["gasch_comparison"]["gasch_three"],
        tfiib_diamide = FIGURES["stress"]["gasch_comparison"]["tfiib_diamide"],
        tfiib_aminoacid = FIGURES["stress"]["gasch_comparison"]["tfiib_aminoacid"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/stress/stress_gasch_comparison.pdf",
    params:
        width = eval(str(FIGURES["stress"]["gasch_comparison"]["width"])),
        height = eval(str(FIGURES["stress"]["gasch_comparison"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/stress_gasch_comparison.R"

rule nitrogen_gene_ontology:
    input:
        nitrogen_ontology = FIGURES["stress"]["nitrogen_gene_ontology"]["data"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/stress/stress_nitrogen_gene_ontology.pdf",
    params:
        width = eval(str(FIGURES["stress"]["nitrogen_gene_ontology"]["width"])),
        height = eval(str(FIGURES["stress"]["nitrogen_gene_ontology"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/stress_nitrogen_gene_ontology.R"

rule tfiib_ridgelines:
    input:
        data_path = FIGURES["stress"]["tfiib_ridgelines"]["data_path"],
        annotations = FIGURES["stress"]["tfiib_ridgelines"]["annotations"],
        diamide_genic = FIGURES["stress"]["tfiib_ridgelines"]["diamide_genic"],
        diamide_intragenic = FIGURES["stress"]["tfiib_ridgelines"]["diamide_intragenic"],
        amino_genic = FIGURES["stress"]["tfiib_ridgelines"]["amino_genic"],
        amino_intragenic = FIGURES["stress"]["tfiib_ridgelines"]["amino_intragenic"],
        nitrogen_genic = FIGURES["stress"]["tfiib_ridgelines"]["nitrogen_genic"],
        nitrogen_intragenic = FIGURES["stress"]["tfiib_ridgelines"]["nitrogen_intragenic"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/stress/stress_tfiib_ridgelines.pdf",
    params:
        width = eval(str(FIGURES["stress"]["tfiib_ridgelines"]["width"])),
        height = eval(str(FIGURES["stress"]["tfiib_ridgelines"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/stress_tfiib_ridgelines.R"

rule tfiib_coverage:
    input:
        data_paths = [v["paths"] for k,v in  FIGURES["stress"]["tfiib_coverage"]["data"].items()],
        transcript_annotation = FIGURES["stress"]["tfiib_coverage"]["transcript_annotation"],
        orf_annotation = FIGURES["stress"]["tfiib_coverage"]["orf_annotation"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/stress/stress_tfiib_coverage.pdf",
    params:
        gene_ids = FIGURES["stress"]["tfiib_coverage"]["data"].keys(),
        upstream = [v["upstream"] for k,v in FIGURES["stress"]["tfiib_coverage"]["data"].items()],
        downstream = [v["downstream"] for k,v in FIGURES["stress"]["tfiib_coverage"]["data"].items()],
        width = eval(str(FIGURES["stress"]["tfiib_coverage"]["width"])),
        height = eval(str(FIGURES["stress"]["tfiib_coverage"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/stress_tfiib_coverage.R"

rule genic_vs_intra:
    input:
        diamide = FIGURES["stress"]["genic_vs_intra"]["diamide"],
        aminoacid = FIGURES["stress"]["genic_vs_intra"]["aminoacid"],
        nitrogen = FIGURES["stress"]["genic_vs_intra"]["nitrogen"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/stress/stress_genic_vs_intra.pdf",
    params:
        width = eval(str(FIGURES["stress"]["genic_vs_intra"]["width"])),
        height = eval(str(FIGURES["stress"]["genic_vs_intra"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/stress_genic_vs_intra.R"

rule promoter_tss_diffexp_summary:
    input:
        genic_path = FIGURES["stress"]["promoter_tss_diffexp_summary"]["genic_path"],
        intra_path = FIGURES["stress"]["promoter_tss_diffexp_summary"]["intra_path"],
        anti_path = FIGURES["stress"]["promoter_tss_diffexp_summary"]["anti_path"],
        inter_path = FIGURES["stress"]["promoter_tss_diffexp_summary"]["inter_path"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/stress/stress_promoter_tss_diffexp_summary.pdf",
    params:
        fdr_cutoff_tss = config["tss_fdr_cutoff"],
        width = eval(str(FIGURES["stress"]["promoter_tss_diffexp_summary"]["width"])),
        height = eval(str(FIGURES["stress"]["promoter_tss_diffexp_summary"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/stress_promoter_tss_diffexp_summary.R"

rule promoter_tss_expression:
    input:
        genic_path = FIGURES["stress"]["promoter_tss_expression"]["genic_path"],
        intra_path = FIGURES["stress"]["promoter_tss_expression"]["intra_path"],
        counts_path = FIGURES["stress"]["promoter_tss_expression"]["counts_path"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/stress/stress_promoter_tss_expression.pdf",
    params:
        width = eval(str(FIGURES["stress"]["promoter_tss_expression"]["width"])),
        height = eval(str(FIGURES["stress"]["promoter_tss_expression"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/stress_promoter_tss_expression.R"

rule promoter_tss_polyenrichment:
    input:
        matched_peaks_genic = FIGURES["stress"]["promoter_tss_polyenrichment"]["matched_peaks_genic"],
        matched_peaks_intra = FIGURES["stress"]["promoter_tss_polyenrichment"]["matched_peaks_intra"],
        polyenrichment_genic = FIGURES["stress"]["promoter_tss_polyenrichment"]["polyenrichment_genic"],
        polyenrichment_intra = FIGURES["stress"]["promoter_tss_polyenrichment"]["polyenrichment_intra"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/stress/stress_promoter_tss_polyenrichment.pdf",
    params:
        fdr_cutoff_tss = config["tss_fdr_cutoff"],
        width = eval(str(FIGURES["stress"]["promoter_tss_polyenrichment"]["width"])),
        height = eval(str(FIGURES["stress"]["promoter_tss_polyenrichment"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/stress_promoter_tss_polyenrichment.R"

rule dsk2_summary:
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
        pdf = "figures/stress/stress_dsk2_summary.pdf",
    params:
        gene_ids = FIGURES["stress"]["dsk2_summary"]["gene_ids"],
        up_distances = FIGURES["stress"]["dsk2_summary"]["up_distances"],
        down_distances = FIGURES["stress"]["dsk2_summary"]["down_distances"],
        modification = FIGURES["stress"]["dsk2_summary"]["modification"],
        refpointlabel = FIGURES["stress"]["dsk2_summary"]["refpointlabel"],
        trim_pct = FIGURES["stress"]["dsk2_summary"]["trim_pct"],
        x_label_dist = FIGURES["stress"]["dsk2_summary"]["x_label_dist"],
        width = eval(str(FIGURES["stress"]["dsk2_summary"]["width"])),
        height = eval(str(FIGURES["stress"]["dsk2_summary"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/stress_dsk2_summary.R"

rule dsk2_pace_northern:
    input:
        dsk_blot = FIGURES["stress"]["dsk2_pace_northern"]["dsk_blot"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/stress/stress_dsk2_pace_northern.pdf",
    params:
        width = eval(str(FIGURES["stress"]["dsk2_pace_northern"]["width"])),
        height = eval(str(FIGURES["stress"]["dsk2_pace_northern"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/stress_dsk2_pace_northern.R"

rule diamide_fitnesscomp:
    input:
        data_path = FIGURES["stress"]["diamide_fitnesscomp"]["data_path"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/stress/stress_diamide_fitnesscomp.pdf",
    params:
        width = eval(str(FIGURES["stress"]["diamide_fitnesscomp"]["width"])),
        height = eval(str(FIGURES["stress"]["diamide_fitnesscomp"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/stress_diamide_fitnesscomp.R"

rule interyeast_intragenic:
    input:
        mikatae = FIGURES["stress"]["interyeast_intragenic"]["mikatae"],
        uvarum = FIGURES["stress"]["interyeast_intragenic"]["uvarum"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/stress/stress_interyeast_intragenic.pdf",
    params:
        width = eval(str(FIGURES["stress"]["interyeast_intragenic"]["width"])),
        height = eval(str(FIGURES["stress"]["interyeast_intragenic"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/stress_interyeast_intragenic.R"

rule dsk2_interyeast:
    input:
        transcripts_scer = FIGURES["stress"]["dsk2_interyeast"]["transcripts_scer"],
        transcripts_smik = FIGURES["stress"]["dsk2_interyeast"]["transcripts_smik"],
        transcripts_sbay = FIGURES["stress"]["dsk2_interyeast"]["transcripts_sbay"],
        orfs_scer = FIGURES["stress"]["dsk2_interyeast"]["orfs_scer"],
        orfs_smik = FIGURES["stress"]["dsk2_interyeast"]["orfs_smik"],
        orfs_sbay = FIGURES["stress"]["dsk2_interyeast"]["orfs_sbay"],
        pace_scer = FIGURES["stress"]["dsk2_interyeast"]["pace_scer"],
        pace_smik = FIGURES["stress"]["dsk2_interyeast"]["pace_smik"],
        pace_sbay = FIGURES["stress"]["dsk2_interyeast"]["pace_sbay"],
        peak_scer = FIGURES["stress"]["dsk2_interyeast"]["peak_scer"],
        peak_smik = FIGURES["stress"]["dsk2_interyeast"]["peak_smik"],
        peak_sbay = FIGURES["stress"]["dsk2_interyeast"]["peak_sbay"],
        tss_scer = FIGURES["stress"]["dsk2_interyeast"]["tss_scer"],
        tss_smik = FIGURES["stress"]["dsk2_interyeast"]["tss_smik"],
        tss_sbay = FIGURES["stress"]["dsk2_interyeast"]["tss_sbay"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/stress/stress_dsk2_interyeast.pdf",
    params:
        width = eval(str(FIGURES["stress"]["dsk2_interyeast"]["width"])),
        height = eval(str(FIGURES["stress"]["dsk2_interyeast"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/stress_dsk2_interyeast.R"

rule dsk2_interyeast_northern:
    input:
        dsk_scer = FIGURES["stress"]["dsk2_interyeast_northern"]["dsk_scer"],
        snr_scer = FIGURES["stress"]["dsk2_interyeast_northern"]["snr_scer"],
        hsp_scer = FIGURES["stress"]["dsk2_interyeast_northern"]["hsp_scer"],
        dsk_smik = FIGURES["stress"]["dsk2_interyeast_northern"]["dsk_smik"],
        snr_smik = FIGURES["stress"]["dsk2_interyeast_northern"]["snr_smik"],
        hsp_smik = FIGURES["stress"]["dsk2_interyeast_northern"]["hsp_smik"],
        dsk_sbay = FIGURES["stress"]["dsk2_interyeast_northern"]["dsk_sbay"],
        snr_sbay = FIGURES["stress"]["dsk2_interyeast_northern"]["snr_sbay"],
        hsp_sbay = FIGURES["stress"]["dsk2_interyeast_northern"]["hsp_sbay"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/stress/stress_dsk2_interyeast_northern.pdf",
    params:
        width = eval(str(FIGURES["stress"]["dsk2_interyeast_northern"]["width"])),
        height = eval(str(FIGURES["stress"]["dsk2_interyeast_northern"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/stress_dsk2_interyeast_northern.R"

