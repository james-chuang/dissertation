
main = function(theme_spec,
                fonts_path,
                heatmap_scripts,
                rnaseq_data, annotation,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)
    source(heatmap_scripts)

    ttf_import(fonts_path)
    loadfonts()

    sample_list = c("non-depleted-1", "non-depleted-2", "depleted-1", "depleted-2")

    plot = plot_heatmap(data_path = rnaseq_data,
                        sample_list = sample_list,
                        anno_path = annotation,
                        cps_dist = 0.3,
                        experiment = "spt5",
                        max_length=2.5,
                        cutoff_pct=0.9,
                        add_ylabel=TRUE,
                        y_label = "nonoverlapping coding genes",
                        colorbar_title="antisense RNA-seq signal")

    ggsave(pdf_out, plot=plot, width=fig_width, height=fig_height, units="in")
    embed_fonts(pdf_out)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path = snakemake@input[["fonts_path"]],
     heatmap_scripts = snakemake@input[["heatmap_scripts"]],
     rnaseq_data = snakemake@input[["rnaseq_data"]],
     annotation = snakemake@input[["annotation"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

