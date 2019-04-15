#!/usr/bin/env python

rule five_netseq_metagenes:
    input:
        netseq_data = FIGURES["five"]["netseq_metagene"]["netseq_data"],
        chipseq_data = FIGURES["five"]["netseq_metagene"]["chipseq_data"],
        annotation = FIGURES["five"]["netseq_metagene"]["annotation"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/five/five_netseq_metagene.pdf",
    params:
        height = eval(str(FIGURES["five"]["netseq_metagene"]["height"])),
        width = eval(str(FIGURES["five"]["netseq_metagene"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_netseq_metagene.R"



