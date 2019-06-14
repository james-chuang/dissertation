#!/usr/bin/env python

configfile: "config.yaml"

FIGURES = config["figures"]

include: "rules/six_figures.smk"
include: "rules/five_figures.smk"
include: "rules/stress_figures.smk"
include: "rules/presentation_figures.smk"

onsuccess:
    shell("(bash compile.sh) > compile.log")

rule target:
    input:
        ".fonts_registered.txt",
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
        "figures/stress/stress_genic_vs_intra.pdf",
        "figures/stress/stress_promoter_tss_diffexp_summary.pdf",
        "figures/stress/stress_promoter_tss_expression.pdf",
        "figures/stress/stress_promoter_tss_polyenrichment.pdf",
        "figures/stress/stress_dsk2_summary.pdf",
        "figures/stress/stress_dsk2_pace_northern.pdf",
        "figures/stress/stress_diamide_fitnesscomp.pdf",
        "figures/stress/stress_interyeast_intragenic.pdf",
        "figures/stress/stress_dsk2_interyeast.pdf",
        "figures/stress/stress_dsk2_interyeast_northern.pdf",
        "figures/presentation/presentation_six_spt6_western.pdf",
        "figures/presentation/presentation_six_gene_diagram.pdf",
        "figures/presentation/presentation_six_aat_assay_comparison.pdf",
        "figures/presentation/presentation_six_tss_seq_heatmaps.pdf",
        "figures/presentation/presentation_six_tss_diffexp_summary.pdf",
        "figures/presentation/presentation_six_tss_expression_levels.pdf",
        "figures/presentation/presentation_six_tfiib_heatmap.pdf",
        expand("figures/presentation/presentation_six_tfiib_spreading_swc{frame}.pdf", frame=["0001", "0002"]),
        "figures/presentation/presentation_six_tss_v_tfiib.pdf",
        "figures/presentation/presentation_six_mnase_heatmaps.pdf",
        expand("figures/presentation/presentation_six_mnase_metagene{frame}.pdf", frame=["0001", "0002"]),
        "figures/presentation/presentation_six_intragenic_mnase_metagenes.pdf",
        "figures/presentation/presentation_six_intragenic_gc.pdf",
        expand("figures/presentation/presentation_six_intragenic_tata{frame}.pdf", frame=["0001", "0002"]),
        "figures/presentation/presentation_five_spt5_depletion.pdf",
        "figures/presentation/presentation_five_netseq_meta.pdf",
        "figures/presentation/presentation_five_rnapii_phosphomark_enrichment.pdf",
        "figures/presentation/presentation_five_rnaseq_metagene.pdf",
        "figures/presentation/presentation_five_rnaseq_heatmaps.pdf",
        "figures/presentation/presentation_five_antisense_heatmaps.pdf",
        "figures/presentation/presentation_five_mnase_metagene.pdf",
        "figures/presentation/presentation_stress_tfiib_ridgelines.pdf",
        "figures/presentation/presentation_stress_tfiib_coverage.pdf",
        "figures/presentation/presentation_stress_promoter_tss_expression.pdf",
        "figures/presentation/presentation_stress_promoter_tss_polyenrichment.pdf",
        "figures/presentation/presentation_stress_dsk2_summary.pdf",
        "figures/presentation/presentation_stress_dsk2_pace_northern.pdf",
        "figures/presentation/presentation_stress_diamide_fitnesscomp.pdf",

rule register_fonts:
    input:
        fonts_path = config["fonts_path"],
    output:
        output_path = ".fonts_registered.txt"
    conda:
        "envs/plot.yaml"
    script:
        "scripts/register_fonts.R"

