
main = function(theme_spec,
                fonts_path,
                heatmap_scripts,
                tfiib_data,
                annotation,
                fig_width,
                fig_height,
                pdf_out){
    source(theme_spec)
    source(heatmap_scripts)

    ttf_import(fonts_path)
    loadfonts()

    sample_list = c("WT-37C-1", "WT-37C-2", "spt6-1004-37C-1", "spt6-1004-37C-2")

    fig_two_a = plot_heatmap(data_path = tfiib_data,
                             sample_list = sample_list,
                             anno_path = annotation,
                             cps_dist = 0.3,
                             max_length=2.5,
                             cutoff_pct=0.85,
                             add_ylabel=TRUE,
                             y_label = "nonoverlapping coding genes",
                             colorbar_title="TFIIB ChIP-nexus protection")

    ggsave(pdf_out, plot=fig_two_a, width=fig_width, height=fig_height, units="in")
    embed_fonts(pdf_out)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path = snakemake@input[["fonts_path"]],
     heatmap_scripts = snakemake@input[["heatmap_scripts"]],
     tfiib_data = snakemake@input[["tfiib_data"]],
     annotation = snakemake@input[["annotation"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

