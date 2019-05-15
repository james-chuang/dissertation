#!/usr/bin/env python

configfile: "config.yaml"

FIGURES = config["figures"]

include: "rules/six_figures.smk"
include: "rules/five_figures.smk"
include: "rules/stress_figures.smk"

onsuccess:
    shell("(bash compile.sh) > compile.log")

rule target:
    input:
        "figures/six/six_spt6_western.pdf",
        "figures/six/six_gene_diagram.pdf",
        "figures/six/six_aat_assay_comparison.pdf",
        "figures/six/six_tss_seq_heatmaps.pdf",
        "figures/six/six_tfiib_heatmap.pdf",
        "figures/six/six_tss_diffexp_summary.pdf",
        "figures/six/six_tss_expression_levels.pdf",
        "figures/six/six_intragenic_genes_bvenn.pdf",
        "figures/six/six_tfiib_spreading_ssa4.pdf",
        "figures/six/six_tss_v_tfiib.pdf",
        "figures/six/six_mnase_metagene.pdf",
        "figures/six/six_global_nuc_occ_fuzz.pdf",
        "figures/six/six_mnase_som.pdf",
        "figures/six/six_mnase_heatmaps.pdf",
        "figures/six/six_intragenic_mnase_metagenes.pdf",
        "figures/six/six_tss_seqlogos.pdf",
        "figures/six/six_intragenic_tata.pdf",
        "figures/six/six_meme_motifs.pdf",
        "figures/five/five_summary_metagenes.pdf",
        "figures/five/five_rnapii_phosphomark_enrichment.pdf",
        "figures/five/five_tss_vs_rna.pdf",
        "figures/five/five_rnaseq_metagene.pdf",
        "figures/five/five_rnaseq_heatmaps.pdf",
        "figures/five/five_antisense_heatmaps.pdf",
        "figures/five/five_tss_diffexp_summary.pdf",
        "figures/five/five_tss_expression_levels.pdf",
        "figures/five/five_meme_motifs.pdf",
        "figures/five/five_mnase_metagene.pdf",
        "figures/five/five_nuc_fuzz.pdf",
        "figures/five/five_antisense_mnase_metagene.pdf",
        "figures/stress/stress_gasch_comparison.pdf",
        "figures/stress/stress_nitrogen_gene_ontology.pdf",
        "figures/stress/stress_tfiib_ridgelines.pdf",
        "figures/stress/stress_tfiib_coverage.pdf",
        "figures/stress/stress_promoter_tss_diffexp_summary.pdf",
        "figures/stress/stress_promoter_tss_expression.pdf",
        "figures/stress/stress_promoter_tss_polyenrichment.pdf",
