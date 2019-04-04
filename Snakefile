#!/usr/bin/env python

configfile: "config.yaml"

FIGURES = config["figures"]

include:
    "rules/stress_figures.smk"

rule target:
    input:
        "figures/stress/stress_tfiib_ridgelines.pdf",

