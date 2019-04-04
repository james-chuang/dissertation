#!/usr/bin/env python

configfile: "config.yaml"

FIGURES = config["figures"]

include:
    "rules/stress_figures.smk"

onsuccess:
    shell("(bash compile.sh) > compile.log")


rule target:
    input:
        "figures/stress/stress_tfiib_ridgelines.pdf",
        "figures/stress/stress_tfiib_coverage.pdf",
        "figures/stress/stress_promoter_tss_diffexp_summary.pdf",
        "figures/stress/stress_promoter_tss_expression.pdf",
        "figures/stress/stress_promoter_tss_polyenrichment.pdf",
