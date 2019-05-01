#!/usr/bin/env python

rule spt6_western:
    input:
        data_path = FIGURES["six"]["spt6_western"]["data_path"],
        blot_path = FIGURES["six"]["spt6_western"]["blot_path"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_spt6_western.pdf",
    params:
        height = eval(str(FIGURES["six"]["spt6_western"]["height"])),
        width = eval(str(FIGURES["six"]["spt6_western"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_spt6_western.R"

rule six_gene_diagram:
    input:
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"]
    output:
        pdf = "figures/six/six_gene_diagram.pdf",
    params:
        width = eval(str(FIGURES["six"]["gene_diagram"]["width"])),
        height = eval(str(FIGURES["six"]["gene_diagram"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_gene_diagram.R"

rule six_aat_assay_comparison:
    input:
        data_path = FIGURES["six"]["aat_assay_comparison"]["data_path"],
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"]
    output:
        pdf = "figures/six/six_aat_assay_comparison.pdf",
    params:
        width = eval(str(FIGURES["six"]["aat_assay_comparison"]["width"])),
        height = eval(str(FIGURES["six"]["aat_assay_comparison"]["height"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_aat2_assay_comparisons.R"

rule six_tss_seq_heatmaps:
    input:
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"],
        heatmap_scripts = "scripts/plot_heatmap.R",
        annotation = FIGURES["six"]["tss_seq_heatmaps"]["annotation"],
        tss_sense = FIGURES["six"]["tss_seq_heatmaps"]["tss_sense"],
        tss_antisense = FIGURES["six"]["tss_seq_heatmaps"]["tss_antisense"]
    output:
        pdf = "figures/six/six_tss_seq_heatmaps.pdf"
    params:
        height = eval(str(FIGURES["six"]["tss_seq_heatmaps"]["height"])),
        width = eval(str(FIGURES["six"]["tss_seq_heatmaps"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_tss_seq_heatmaps.R"

rule six_tfiib_heatmap:
    input:
        fonts_path = config["fonts_path"],
        theme = config["theme_spec"],
        heatmap_scripts = "scripts/plot_heatmap.R",
        annotation = FIGURES["six"]["tfiib_heatmap"]["annotation"],
        tfiib_data = FIGURES["six"]["tfiib_heatmap"]["tfiib_data"],
    output:
        pdf = "figures/six/six_tfiib_heatmap.pdf"
    params:
        height = eval(str(FIGURES["six"]["tfiib_heatmap"]["height"])),
        width = eval(str(FIGURES["six"]["tfiib_heatmap"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_tfiib_heatmap.R"

rule six_tss_diffexp_summary:
    input:
        in_genic = FIGURES["six"]["tss_diffexp_summary"]["genic"],
        in_intra = FIGURES["six"]["tss_diffexp_summary"]["intragenic"],
        in_anti = FIGURES["six"]["tss_diffexp_summary"]["antisense"],
        in_inter = FIGURES["six"]["tss_diffexp_summary"]["intergenic"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_tss_diffexp_summary.pdf",
    params:
        fdr_cutoff_tss = config["tss_fdr_cutoff"],
        height = eval(str(FIGURES["six"]["tss_diffexp_summary"]["height"])),
        width = eval(str(FIGURES["six"]["tss_diffexp_summary"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_tss_diffexp_summary.R"

rule six_tss_expression_levels:
    input:
        tss_genic = FIGURES["six"]["tss_expression_levels"]["tss_genic"],
        tss_intragenic = FIGURES["six"]["tss_expression_levels"]["tss_intragenic"],
        tss_antisense = FIGURES["six"]["tss_expression_levels"]["tss_antisense"],
        tss_intergenic = FIGURES["six"]["tss_expression_levels"]["tss_intergenic"],
        # tfiib_genic = FIGURES["figure_one"]["one_d"]["tfiib_genic"],
        # tfiib_intragenic = FIGURES["figure_one"]["one_d"]["tfiib_intragenic"],
        # tfiib_intergenic = FIGURES["figure_one"]["one_d"]["tfiib_intergenic"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_tss_expression_levels.pdf",
    params:
        height = eval(str(FIGURES["six"]["tss_expression_levels"]["height"])),
        width = eval(str(FIGURES["six"]["tss_expression_levels"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_tss_expression_levels.R"

rule intragenic_genes_bvenn:
    input:
        common_names = FIGURES["six"]["intragenic_genes_bvenn"]["common_names"],
        cheung_data = FIGURES["six"]["intragenic_genes_bvenn"]["cheung_data"],
        uwimana_data = FIGURES["six"]["intragenic_genes_bvenn"]["uwimana_data"],
        tss_data = FIGURES["six"]["intragenic_genes_bvenn"]["tss_data"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_intragenic_genes_bvenn.pdf",
    params:
        height = eval(str(FIGURES["six"]["intragenic_genes_bvenn"]["height"])),
        width = eval(str(FIGURES["six"]["intragenic_genes_bvenn"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_intragenic_genes_bvenn.R"

rule six_tfiib_spreading_ssa4:
    input:
        tfiib_data = FIGURES["six"]["tfiib_spreading_ssa4"]["tfiib_data"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_tfiib_spreading_ssa4.pdf",
    params:
        height = eval(str(FIGURES["six"]["tfiib_spreading_ssa4"]["height"])),
        width = eval(str(FIGURES["six"]["tfiib_spreading_ssa4"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_tfiib_spreading_ssa4.R"

rule six_tss_vs_tfiib:
    input:
        genic = FIGURES["six"]["tss_vs_tfiib"]["genic"],
        intragenic = FIGURES["six"]["tss_vs_tfiib"]["intragenic"],
        antisense = FIGURES["six"]["tss_vs_tfiib"]["antisense"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_tss_v_tfiib.pdf",
    params:
        height = eval(str(FIGURES["six"]["tss_vs_tfiib"]["height"])),
        width = eval(str(FIGURES["six"]["tss_vs_tfiib"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_tss_vs_tfiib.R"

rule six_mnase_metagene:
    input:
        mnase_data = FIGURES["six"]["mnase_metagene"]["mnase_data"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_mnase_metagene.pdf",
    params:
        height = eval(str(FIGURES["six"]["mnase_metagene"]["height"])),
        width = eval(str(FIGURES["six"]["mnase_metagene"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_mnase_metagene.R"

rule six_global_nuc_occ_fuzz:
    input:
        wt_mnase_quant = FIGURES["six"]["global_nuc_occ_fuzz"]["wt_mnase_quant"],
        spt6_mnase_quant = FIGURES["six"]["global_nuc_occ_fuzz"]["spt6_mnase_quant"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_global_nuc_occ_fuzz.pdf",
    params:
        height = eval(str(FIGURES["six"]["global_nuc_occ_fuzz"]["height"])),
        width = eval(str(FIGURES["six"]["global_nuc_occ_fuzz"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_global_nuc_occ_fuzz.R"

rule six_mnase_heatmaps:
    input:
        netseq_data = FIGURES["six"]["mnase_heatmaps"]["netseq_data"],
        mnase_data = FIGURES["six"]["mnase_heatmaps"]["mnase_data"],
        quant_data = FIGURES["six"]["mnase_heatmaps"]["quant_data"],
        annotation = FIGURES["six"]["mnase_heatmaps"]["annotation"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_mnase_heatmaps.pdf",
    params:
        height = eval(str(FIGURES["six"]["mnase_heatmaps"]["height"])),
        width = eval(str(FIGURES["six"]["mnase_heatmaps"]["width"])),
        assay = "NET-seq"
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_mnase_heatmaps.R"

rule six_mnase_som:
    input:
        som_data = FIGURES["six"]["mnase_som"]["som_data"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_mnase_som.pdf",
    params:
        height = eval(str(FIGURES["six"]["mnase_som"]["height"])),
        width = eval(str(FIGURES["six"]["mnase_som"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_mnase_som.R"

rule six_intragenic_mnase_metagenes:
    input:
        mnase_data = FIGURES["six"]["intragenic_mnase_metagenes"]["mnase_data"],
        gc_data = FIGURES["six"]["intragenic_mnase_metagenes"]["gc_data"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_intragenic_mnase_metagenes.pdf",
    params:
        height = eval(str(FIGURES["six"]["intragenic_mnase_metagenes"]["height"])),
        width = eval(str(FIGURES["six"]["intragenic_mnase_metagenes"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_intragenic_mnase_metagenes.R"

rule six_tss_seqlogos:
    input:
        data_paths = FIGURES["six"]["tss_seqlogos"]["data_paths"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_tss_seqlogos.pdf",
    params:
        height = eval(str(FIGURES["six"]["tss_seqlogos"]["height"])),
        width = eval(str(FIGURES["six"]["tss_seqlogos"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_tss_seqlogos.R"

rule six_intragenic_tata:
    input:
        tata_genic_path = FIGURES["six"]["intragenic_tata"]["tata_genic"],
        tata_intra_path = FIGURES["six"]["intragenic_tata"]["tata_intragenic"],
        tata_anti_path = FIGURES["six"]["intragenic_tata"]["tata_antisense"],
        # tata_random_path = FIGURES["six"]["intragenic_tata"]["tata_random"],
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_intragenic_tata.pdf",
    params:
        height = eval(str(FIGURES["six"]["intragenic_tata"]["height"])),
        width = eval(str(FIGURES["six"]["intragenic_tata"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_intragenic_tata.R"

rule six_build_motif_input:
    input:
        intragenic = FIGURES["six"]["meme_motifs"]["intragenic"],
        antisense = FIGURES["six"]["meme_motifs"]["antisense"],
    output:
        intragenic = "figures/six/six_intragenic_motif.meme",
        antisense = "figures/six/six_antisense_motif.meme"
    conda:
        "../envs/meme.yaml"
    shell: """
        meme2meme <(meme-get-motif -a -id GAHRATGAAGAWGADGAHGAT-MEME-1 {input.intragenic}) > {output.intragenic}
        meme2meme <(meme-get-motif -a -id TCWTCDTCNTCDWCTTCTTCWTCDTYWTC-MEME-1 {input.antisense}) \\
                <(meme-get-motif -a -id TATWTAKATATATATATAYWT-MEME-2 {input.antisense}) \\
                <(meme-get-motif -a -id RAARAAADWRAAAAARAARRA-MEME-3 {input.antisense}) > {output.antisense}
        """

rule six_meme_motifs:
    input:
        intragenic = "figures/six/six_intragenic_motif.meme",
        antisense = "figures/six/six_antisense_motif.meme",
        theme = config["theme_spec"],
        fonts_path = config["fonts_path"],
    output:
        pdf = "figures/six/six_meme_motifs.pdf",
    params:
        height = eval(str(FIGURES["six"]["meme_motifs"]["height"])),
        width = eval(str(FIGURES["six"]["meme_motifs"]["width"])),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/six_meme_motifs.R"

