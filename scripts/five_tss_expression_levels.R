
import = function(df, path, category){
    df = read_tsv(path) %>%
        filter(! is.na(log10_padj)) %>%
        mutate(category=category) %>%
        bind_rows(df, .) %>%
        return()
}

main = function(theme_spec,
                fonts_path,
                tss_genic, tss_intragenic, tss_antisense, tss_intergenic,
                # tfiib_genic, tfiib_intragenic, tfiib_intergenic,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)

    ttf_import(fonts_path)
    loadfonts()

    tss_df = tibble() %>%
        import(tss_genic, "genic") %>%
        import(tss_intragenic, "intragenic") %>%
        import(tss_antisense, "antisense") %>%
        import(tss_intergenic, "intergenic") %>%
        mutate(category=fct_inorder(category, ordered=TRUE)) %>%
        select(name, condition_expr, control_expr, category) %>%
        gather(key=condition, value=expression, -c(name, category)) %>%
        mutate(condition = ordered(condition,
                                   levels = c("control_expr", "condition_expr"),
                                   labels = c("non-depleted", "depleted")))

    tss_plot = ggplot(data = tss_df,
                      aes(x=category, y=expression+1,
                          group=interaction(condition, category))) +
        geom_violin(aes(fill=condition),
                    bw = .08,
                    width=1.2,
                    position=position_dodge(width=0.75),
                    size=0.2) +
        geom_boxplot(position=position_dodge(width=0.75),
                     width=0.12,
                     notch=TRUE,
                     outlier.size=0,
                     outlier.stroke=0,
                     size=0.2) +
        scale_x_discrete(expand = c(0,0)) +
                         # limits = c("genic", "intragenic", "antisense", "intergenic", "")) +
        scale_y_log10(name = "normalized counts",
                      breaks = c(10, 1000), labels = c(bquote(10^1), bquote(10^3))) +
        scale_fill_few(guide=guide_legend(direction="vertical",
                                          label.position="left",
                                          # label.hjust=1,
                                          keyheight=unit(6, "pt"))) +
        # ggtitle("expression level of TSS-seq peaks") +
        theme_default +
        theme(axis.title.x = element_blank(),
              panel.grid.major.x = element_blank(),
              legend.position = "top",
              legend.text = element_text(margin = margin(0,0,0,0,"pt")),
              legend.justification = c(0.5, 0.5),
              legend.box.margin = margin(0, 0, -10, 0, "pt"),
              plot.margin = margin(0,0,0,0,"pt"))

    ggsave(pdf_out,
           plot=tss_plot,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path= snakemake@input[["fonts_path"]],
     tss_genic = snakemake@input[["tss_genic"]],
     tss_intragenic = snakemake@input[["tss_intragenic"]],
     tss_antisense = snakemake@input[["tss_antisense"]],
     tss_intergenic = snakemake@input[["tss_intergenic"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

