filter_upregulated_promoters = function(tfiib_match_path, fdr_cutoff_tss, polyenrichment_path, category_id, feature_name_col){
    feature_name_col = enquo(feature_name_col)

    peak_list = read_tsv(tfiib_match_path) %>%
        filter(tss_FDR > -log10(fdr_cutoff_tss) &
                   tss_lfc >= 0 &
                   tfiib_start >= 0) %>%
        distinct(tss_name) %>%
        pull(tss_name)

    read_tsv(polyenrichment_path) %>%
        filter(name %in% peak_list) %>%
        select(1:11,
               feature_name = !!feature_name_col,
               log2_polyenrichment, lfc_SE_polyenrichment, log10_padj_polyenrichment,
               condition_polyenrichment, condition_polyenrichment_SE,
               control_polyenrichment, control_polyenrichment_SE) %>%
        mutate(category = category_id) %>%
        return()
}

main = function(theme_spec,
                fdr_cutoff_tss,
                matched_peaks_genic,
                matched_peaks_intra,
                polyenrichment_genic,
                polyenrichment_intra,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)

    df = filter_upregulated_promoters(matched_peaks_genic,
                                      fdr_cutoff_tss,
                                      polyenrichment_genic,
                                      "genic",
                                      genic_name) %>%
        bind_rows(filter_upregulated_promoters(matched_peaks_intra,
                                      fdr_cutoff_tss,
                                      polyenrichment_intra,
                                      "intragenic",
                                      orf_name)) %>%
        mutate(category = fct_inorder(category)) %>%
        group_by(category) %>%
        mutate(rank = rank(condition_polyenrichment),
               y = rank/max(rank))

    figure_3b = ggplot(data = df) +
        geom_segment(data = df,
                     aes(x=condition_polyenrichment-condition_polyenrichment_SE,
                         xend=condition_polyenrichment+condition_polyenrichment_SE,
                         y=y,
                         yend=y,
                         color=category),
                        alpha=0.4,
                        size=0.7) +
        geom_point(data = df,
                   aes(x=condition_polyenrichment,
                       y=y,
                       color=category),
                   shape=16,
                   alpha=0.95,
                   size=1.2) +
        geom_label(data = df %>%
                       filter(category=="intragenic"),
                   aes(x=condition_polyenrichment-condition_polyenrichment_SE-0.1,
                       y=y,
                       label=feature_name),
                  hjust=1,
                  size=8/72*25.4,
                  label.r = unit(0, "pt"),
                  label.padding = unit(1, "pt"),
                  label.size = NA,
                  family="FreeSans",
                  fontface="italic") +
        scale_y_continuous(name="percentile",
                           expand = c(0,0.03),
                           breaks = scales::pretty_breaks(n=3),
                           labels = function(x) 100*x) +
        scale_x_continuous(name=expression("relative polysome enrichment" %==%
                                               "log"[2] ~ textstyle(frac("polysome RNA",
                                                                       "total RNA"))),
                           expand=c(.06,0)) +
        ggtitle("polysome enrichment in oxidative stress",
                subtitle = "oxidative stress-induced TSSs") +
        scale_color_tableau() +
        theme_default +
        theme(legend.position = c(0.4, 0.7),
              plot.subtitle = element_text(size=8,
                                           margin=margin(t=0, b=-2, unit="pt")),
              plot.margin=margin(0.5,0,0,0,"pt"))

    ggplot2::ggsave(pdf_out,
                    plot=figure_3b,
                    width=fig_width,
                    height=fig_height,
                    units="in",
                    device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     fdr_cutoff_tss = snakemake@params[["fdr_cutoff_tss"]],
     matched_peaks_genic = snakemake@input[["matched_peaks_genic"]],
     matched_peaks_intra = snakemake@input[["matched_peaks_intra"]],
     polyenrichment_genic = snakemake@input[["polyenrichment_genic"]],
     polyenrichment_intra = snakemake@input[["polyenrichment_intra"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

