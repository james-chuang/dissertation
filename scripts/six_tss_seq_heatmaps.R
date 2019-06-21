
main = function(theme_spec,
                heatmap_scripts,
                annotation,
                sense_tss_data,
                antisense_tss_data,
                fig_width,
                fig_height,
                pdf_out){
    source(theme_spec)
    source(heatmap_scripts)

    sample_list = c("WT-37C-1", "WT-37C-2", "spt6-1004-37C-1", "spt6-1004-37C-2")
    cps_dist = 0.3
    max_length = 2.5
    cutoff_pct = 0.92

    sense_heatmap = plot_heatmap(data_path=sense_tss_data,
                                 sample_list=sample_list,
                                 anno_path=annotation,
                                 cps_dist=cps_dist,
                                 max_length=max_length,
                                 cutoff_pct=cutoff_pct,
                                 add_ylabel=TRUE,
                                 y_label="nonoverlapping coding genes",
                                 palette="magma",
                                 colorbar_title="sense TSS-seq signal")

    anti_heatmap = plot_heatmap(data_path=antisense_tss_data,
                                sample_list=sample_list,
                                anno_path=annotation,
                                cps_dist=cps_dist,
                                max_length=max_length,
                                cutoff_pct=cutoff_pct,
                                add_ylabel=FALSE,
                                palette="magma",
                                colorbar_title="antisense TSS-seq signal")
    figure_two = arrangeGrob(sense_heatmap, anti_heatmap, nrow=1)

    ggsave(pdf_out,
           plot=figure_two,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     heatmap_scripts = snakemake@input[["heatmap_scripts"]],
     annotation = snakemake@input[["annotation"]],
     sense_tss_data = snakemake@input[["tss_sense"]],
     antisense_tss_data = snakemake@input[["tss_antisense"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

