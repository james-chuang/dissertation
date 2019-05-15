#!/usr/bin/env python

rule gasch_comparison:
    input:
        gasch_one = FIGURES["stress"]["gasch_comparison"]["gasch_one"],
        gasch_three = FIGURES["stress"]["gasch_comparison"]["gasch_three"],
        tfiib_diamide = FIGURES["stress"]["gasch_comparison"]["tfiib_diamide"],
        tfiib_aminoacid = FIGURES["stress"]["gasch_comparison"]["tfiib_aminoacid"],
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"]
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
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"]
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
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"]
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
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"]
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

rule promoter_tss_diffexp_summary:
    input:
        genic_path = FIGURES["stress"]["promoter_tss_diffexp_summary"]["genic_path"],
        intra_path = FIGURES["stress"]["promoter_tss_diffexp_summary"]["intra_path"],
        anti_path = FIGURES["stress"]["promoter_tss_diffexp_summary"]["anti_path"],
        inter_path = FIGURES["stress"]["promoter_tss_diffexp_summary"]["inter_path"],
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"]
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
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"]
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
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"]
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


