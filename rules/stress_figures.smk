#!/usr/bin/env python

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


