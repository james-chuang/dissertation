
main = function(theme_spec,
                heatmap_scripts,
                rnaseq_data, annotation,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)
    source(heatmap_scripts)

    sample_list = c("non-depleted-1", "non-depleted-2", "depleted-1", "depleted-2")

    plot = plot_heatmap(data_path = rnaseq_data,
                        sample_list = sample_list,
                        anno_path = annotation,
                        cps_dist = 0.3,
                        experiment = "spt5",
                        max_length=2.5,
                        cutoff_pct=0.93,
                        add_ylabel=TRUE,
                        y_label = "nonoverlapping coding genes",
                        palette="magma",
                        colorbar_title="antisense RNA-seq",
                        theme="presentation")

    ggsave(pdf_out,
           plot=plot,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     heatmap_scripts = snakemake@input[["heatmap_scripts"]],
     rnaseq_data = snakemake@input[["rnaseq_data"]],
     annotation = snakemake@input[["annotation"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

