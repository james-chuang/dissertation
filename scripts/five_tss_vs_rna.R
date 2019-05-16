
main = function(theme_spec = "thesis_theme.R",
                fonts_path,
                tss_data = "depleted-v-non-depleted_tss-seq-spikenorm-peaks-diffexp-results-genic-all.tsv",
                rna_data = "depleted-v-non-depleted_rnaseq-spikenorm-transcripts-diffexp-results-genic-all.tsv",
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)

    ttf_import(fonts_path)
    loadfonts()

    df = read_tsv(tss_data) %>%
        inner_join(read_tsv(rna_data),
                   by = c("genic_name"="name"),
                   suffix = c("_tss","_rna"))

    df %>%
        summarise(tss_median=median(log2_foldchange_tss),
                  rna_median=median(log2_foldchange_rna)) %>%
        print()

    plot = ggplot(data = df) +
        geom_hline(yintercept=0,
                   size=0.5,
                   color="grey70") +
        geom_vline(xintercept=0,
                   size=0.5,
                   color="grey70") +
        geom_abline(slope=1,
                    intercept=0,
                    size=0.5,
                    color="grey70") +
        annotate(geom="label",
                 x = min(df[["log2_foldchange_tss"]]),
                 y = max(df[["log2_foldchange_rna"]]),
                 label = paste("n =", nrow(df)),
                 size = 8/72*25.4,
                 hjust = 0,
                 vjust = 1,
                 label.r = unit(0, "pt"),
                 label.padding = unit(2, "pt"),
                 label.size = NA,
                 family = "FreeSans") +
        stat_bin_hex(geom="point",
                     aes(x=log2_foldchange_tss,
                         y=log2_foldchange_rna,
                         color=..count..),
                     binwidth=c(0.1,0.1),
                     alpha=0.7,
                     shape=16,
                     size=0.5) +
        scale_x_continuous(name = expression("TSS-seq" ~ log[2] ~ textstyle(frac("depleted", "non-depleted"))),
                           breaks = scales::pretty_breaks(n=6)) +
        scale_y_continuous(name = expression(atop("RNA-seq", log[2] ~ textstyle(frac("depleted", "non-depleted")))),
                           breaks = scales::pretty_breaks(n=5)) +
        coord_equal() +
        scale_color_viridis(option="inferno",
                            guide=FALSE) +
        scale_fill_viridis(option="inferno",
                            guide=FALSE) +
        ggtitle("genic transcripts") +
        theme_default +
        theme(axis.title.y = element_text(angle = 0,
                                          vjust = 0.5),
              plot.margin = margin(r=10,
                                   unit="pt"))

    ggsave(pdf_out,
           plot=plot,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path=snakemake@input[["fonts_path"]],
     tss_data = snakemake@input[["tss_data"]],
     rna_data = snakemake@input[["rna_data"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

