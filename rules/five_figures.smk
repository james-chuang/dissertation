#!/usr/bin/env python

rule five_summary_metagenes:
    input:
        netseq_data = FIGURES["five"]["summary_metagenes"]["netseq_data"],
        rnapii_data = FIGURES["five"]["summary_metagenes"]["rnapii_data"],
        spt5_data = FIGURES["five"]["summary_metagenes"]["spt5_data"],
        annotation = FIGURES["five"]["summary_metagenes"]["annotation"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/five/five_summary_metagenes.pdf",
    params:
        height = eval(str(FIGURES["five"]["summary_metagenes"]["height"])),
        width = eval(str(FIGURES["five"]["summary_metagenes"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_summary_metagenes.R"

rule five_rnapii_phosphomark_enrichment:
    input:
        data_paths = FIGURES["five"]["rnapii_phosphomark_enrichment"]["data_paths"],
        annotation = FIGURES["five"]["rnapii_phosphomark_enrichment"]["annotation"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/five/five_rnapii_phosphomark_enrichment.pdf",
    params:
        height = eval(str(FIGURES["five"]["rnapii_phosphomark_enrichment"]["height"])),
        width = eval(str(FIGURES["five"]["rnapii_phosphomark_enrichment"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_rnapii_phosphomark_enrichment.R"

rule five_tss_vs_rna:
    input:
        tss_data = FIGURES["five"]["tss_vs_rna"]["tss_data"],
        rna_data = FIGURES["five"]["tss_vs_rna"]["rna_data"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/five/five_tss_vs_rna.pdf",
    params:
        height = eval(str(FIGURES["five"]["tss_vs_rna"]["height"])),
        width = eval(str(FIGURES["five"]["tss_vs_rna"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_tss_vs_rna.R"

rule five_rnaseq_metagene:
    input:
        data_path = FIGURES["five"]["rnaseq_metagene"]["data_path"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/five/five_rnaseq_metagene.pdf",
    params:
        height = eval(str(FIGURES["five"]["rnaseq_metagene"]["height"])),
        width = eval(str(FIGURES["five"]["rnaseq_metagene"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_rnaseq_metagene.R"


rule five_rnaseq_heatmaps:
    input:
        rnaseq_data = FIGURES["five"]["rnaseq_heatmaps"]["rnaseq_data"],
        annotation = FIGURES["five"]["rnaseq_heatmaps"]["annotation"],
        theme = config["theme_spec"],
        heatmap_scripts = "scripts/plot_heatmap.R",
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/five/five_rnaseq_heatmaps.pdf",
    params:
        height = eval(str(FIGURES["five"]["rnaseq_heatmaps"]["height"])),
        width = eval(str(FIGURES["five"]["rnaseq_heatmaps"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_rnaseq_heatmaps.R"

rule five_antisense_heatmaps:
    input:
        theme = config["theme_spec"],
        tssseq_data = FIGURES["five"]["antisense_heatmaps"]["tssseq_data"],
        rnaseq_data = FIGURES["five"]["antisense_heatmaps"]["rnaseq_data"],
        netseq_data = FIGURES["five"]["antisense_heatmaps"]["netseq_data"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/five/five_antisense_heatmaps.pdf"
    params:
        height = eval(str(FIGURES["five"]["antisense_heatmaps"]["height"])),
        width = eval(str(FIGURES["five"]["antisense_heatmaps"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_antisense_heatmaps.R"

rule five_tss_diffexp_summary:
    input:
        in_genic = FIGURES["five"]["tss_diffexp_summary"]["genic"],
        in_intra = FIGURES["five"]["tss_diffexp_summary"]["intragenic"],
        in_anti = FIGURES["five"]["tss_diffexp_summary"]["antisense"],
        in_inter = FIGURES["five"]["tss_diffexp_summary"]["intergenic"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/five/five_tss_diffexp_summary.pdf",
    params:
        fdr_cutoff_tss = config["tss_fdr_cutoff"],
        height = eval(str(FIGURES["five"]["tss_diffexp_summary"]["height"])),
        width = eval(str(FIGURES["five"]["tss_diffexp_summary"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_tss_diffexp_summary.R"

rule five_tss_expression_levels:
    input:
        tss_genic = FIGURES["five"]["tss_expression_levels"]["tss_genic"],
        tss_intragenic = FIGURES["five"]["tss_expression_levels"]["tss_intragenic"],
        tss_antisense = FIGURES["five"]["tss_expression_levels"]["tss_antisense"],
        tss_intergenic = FIGURES["five"]["tss_expression_levels"]["tss_intergenic"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/five/five_tss_expression_levels.pdf",
    params:
        height = eval(str(FIGURES["five"]["tss_expression_levels"]["height"])),
        width = eval(str(FIGURES["five"]["tss_expression_levels"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_tss_expression_levels.R"

rule five_build_motif_input:
    input:
        antisense = FIGURES["five"]["meme_motifs"]["antisense"],
    output:
        antisense = "figures/five/five_antisense_motif.meme"
    conda:
        "../envs/meme.yaml"
    shell: """
        meme2meme <(meme-get-motif -a -id HTCTTCYTCWTCYWCHTCHTC-MEME-1 {input.antisense}) \\
                <(meme-get-motif -a -id RDRWARWNRAAARAAAAGAAARAAAADM-MEME-2 {input.antisense}) \\
                <(meme-get-motif -a -id CASCAVCVSNWNCWMCDHCAVYDBCWYCRH-MEME-3 {input.antisense}) > {output.antisense}
        """

rule five_meme_motifs:
    input:
        antisense = "figures/five/five_antisense_motif.meme",
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/five/five_meme_motifs.pdf",
    params:
        height = eval(str(FIGURES["five"]["meme_motifs"]["height"])),
        width = eval(str(FIGURES["five"]["meme_motifs"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_meme_motifs.R"

rule five_mnase_metagene:
    input:
        theme = config["theme_spec"],
        mnase_data = FIGURES["five"]["mnase_metagene"]["mnase_data"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/five/five_mnase_metagene.pdf"
    params:
        height = eval(str(FIGURES["five"]["mnase_metagene"]["height"])),
        width = eval(str(FIGURES["five"]["mnase_metagene"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_mnase_metagene.R"

rule five_antisense_mnase_metagene:
    input:
        mnase_data = FIGURES["five"]["antisense_mnase_metagene"]["mnase_data"],
        gc_data = FIGURES["five"]["antisense_mnase_metagene"]["gc_data"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/five/five_antisense_mnase_metagene.pdf",
    params:
        height = eval(str(FIGURES["five"]["antisense_mnase_metagene"]["height"])),
        width = eval(str(FIGURES["five"]["antisense_mnase_metagene"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_antisense_mnase_metagene.R"

rule five_nuc_fuzz:
    input:
        wt_mnase_quant = FIGURES["five"]["nuc_fuzz"]["wt_mnase_quant"],
        mut_mnase_quant = FIGURES["five"]["nuc_fuzz"]["mut_mnase_quant"],
        theme = config["theme_spec"],
        fonts = ".fonts_registered.txt",
    output:
        pdf = "figures/five/five_nuc_fuzz.pdf",
    params:
        height = eval(str(FIGURES["five"]["nuc_fuzz"]["height"])),
        width = eval(str(FIGURES["five"]["nuc_fuzz"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/five_nuc_fuzz.R"


