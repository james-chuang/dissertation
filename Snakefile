#!/usr/bin/env python

configfile: "config.yaml"

FIGURES = config["figures"]

include: "rules/stress_figures.smk"
include: "rules/six_figures.smk"

onsuccess:
    shell("(bash compile.sh) > compile.log")

rule target:
    input:
        "figures/six/six_gene_diagram.pdf",
        "figures/six/six_aat_assay_comparison.pdf",
        "figures/six/six_tss_seq_heatmaps.pdf",
        "figures/six/six_tfiib_heatmap.pdf",
        "figures/six/six_tss_diffexp_summary.pdf",
        "figures/six/six_tss_expression_levels.pdf",
        "figures/stress/stress_tfiib_ridgelines.pdf",
        "figures/stress/stress_tfiib_coverage.pdf",
        "figures/stress/stress_promoter_tss_diffexp_summary.pdf",
        "figures/stress/stress_promoter_tss_expression.pdf",
        "figures/stress/stress_promoter_tss_polyenrichment.pdf",
