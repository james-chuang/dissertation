import = function(path){
    read_tsv(path,
             col_names = c("group", "sample", "annotation",
                           "assay", "index", "position", "signal")) %>%
        group_by(group, annotation, assay, index, position) %>%
        summarise(signal = mean(signal)) %>%
        group_by(group, annotation, assay) %>%
        complete(index, position, fill=list(signal=0)) %>%
        ungroup() %>%
        mutate(group = ordered(group,
                               levels = c("non-depleted",
                                          "depleted"))) %>%
        return()
}

plot_heatmap = function(df, cutoff, colorbar_title, y_title){
    ggplot(data = df, aes(x=position, y=index, fill=signal)) +
        geom_raster(interpolate=FALSE) +
        facet_grid(.~group,
                   switch="y") +
        scale_fill_viridis(option="magma",
                           limits = c(NA, quantile(df[["signal"]], cutoff)),
                           oob = scales::squish,
                           breaks = scales::pretty_breaks(n=2),
                           name = colorbar_title,
                           guide = guide_colorbar(title.position = "top",
                                                  barwidth = 8,
                                                  barheight = 0.3,
                                                  title.hjust=0.5)) +
        scale_x_continuous(breaks = scales::pretty_breaks(n=2),
                           labels = function(x){case_when(near(x, 0) ~ "genic\nTSS",
                                                          near(x, 1) ~ paste(x, "kb"),
                                                          TRUE ~ as.character(x))},
                           expand = c(0, 0)) +
        scale_y_reverse(breaks = function(x){seq(500, max(x), 500)},
                        expand = c(0,0),
                        name = y_title) +
        theme_heatmap_presentation +
        theme(strip.text = element_text(),
              strip.text.x = element_text(vjust=0,
                                          margin = margin(t=-2,
                                                          b=1,
                                                          unit="pt")),
              axis.ticks = element_line(),
              strip.placement = "outside",
              axis.ticks.length = unit(2, "pt"),
              plot.margin = margin(0, 0, 0, 0, "pt"))
}

main = function(theme_spec,
                tss_path,
                # netseq_path,
                rnaseq_path,
                fig_width,
                fig_height,
                pdf_out){
    source(theme_spec)
    library(cowplot)

    tss_df = import(tss_path) %>%
        filter(position <= 1.5)
    rnaseq_df = import(rnaseq_path) %>%
        filter(position <= 1.5)
    # netseq_df = import(netseq_path)

    tss_plot = plot_heatmap(df=tss_df, cutoff=0.995,
                 colorbar_title="antisense TSS-seq",
                 y_title = "932 upregulated antisense TSSs")
    rnaseq_plot = plot_heatmap(df=rnaseq_df, cutoff=0.98,
                 colorbar_title="antisense RNA-seq",
                 y_title = "")
    # netseq_plot = plot_heatmap(df=netseq_df, cutoff=0.96,
    #              colorbar_title="antisense NET-seq",
    #              y_title = "")

    plot = plot_grid(tss_plot,
                     rnaseq_plot,
                     # netseq_plot,
                     align="hv", axis="tbrl", nrow=1)

    ggplot2::ggsave(pdf_out,
                    plot=plot,
                    width=fig_width,
                    height=fig_height,
                    units="cm",
                    device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     tss_path = snakemake@input[["tssseq_data"]],
     # netseq_path = snakemake@input[["netseq_data"]],
     rnaseq_path = snakemake@input[["rnaseq_data"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

