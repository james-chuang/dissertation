
main = function(theme_spec,
                diamide_path,
                aminoacid_path,
                nitrogen_path,
                fig_width, fig_height,
                svg_out, pdf_out, png_out, grob_out){
    source(theme_spec)

    df = read_tsv(diamide_path) %>%
        mutate(condition="oxidative stress") %>%
        bind_rows(read_tsv(aminoacid_path) %>%
                      mutate(condition="amino acid stress")) %>%
        bind_rows(read_tsv(nitrogen_path) %>%
                      mutate(condition="nitrogen stress")) %>%
        mutate(condition = fct_inorder(condition))

    supp_1c = ggplot(data = df) +
        geom_hline(yintercept = 0,
                   color="grey70") +
        geom_vline(xintercept = 0,
                   color="grey70") +
        geom_segment(aes(x=log2_FC_intragenic-lfc_SE_intragenic,
                         xend=log2_FC_intragenic+lfc_SE_intragenic,
                         y=log2_FC_genic,
                         yend=log2_FC_genic,
                         color=condition),
                     size=0.4,
                     alpha=0.3) +
        geom_segment(aes(y=log2_FC_genic-lfc_SE_genic,
                         yend=log2_FC_genic+lfc_SE_genic,
                         x=log2_FC_intragenic,
                         xend=log2_FC_intragenic,
                         color=condition),
                     size=0.4,
                     alpha=0.3) +
        geom_point(aes(x=log2_FC_intragenic,
                       y=log2_FC_genic,
                       color=condition),
                   shape=16,
                   size=1.2,
                   alpha=0.5) +
        annotate(geom="label",
                 x=-4,
                 y=6,
                 label=cor(df[["log2_FC_genic"]],
                           df[["log2_FC_intragenic"]]) %>%
                     signif(2) %>%
                     paste0("R=", .),
                 label.size=NA,
                 label.r = unit(0, "pt"),
                 hjust=0,
                 alpha=0.9,
                 size=8/72*25.4,
                 family="FreeSans") +
        scale_x_continuous(name=expression("intragenic log"[2] ~ textstyle(frac("condition", "control")))) +
        scale_y_continuous(name=expression(atop("genic", "log"[2] ~ textstyle(frac("condition", "control"))))) +
        scale_color_fivethirtyeight() +
        ggtitle("TFIIB ChIP-nexus") +
        theme_default +
        theme(legend.position = c(0.99, 0.99),
              legend.justification = c(1, 1),
              axis.title.y = element_text(angle=0, vjust=0.5, hjust=1),
              plot.margin = margin(r=4, unit="pt"))

    ggplot2::ggsave(pdf_out,
                    plot=supp_1c,
                    width=fig_width,
                    height=fig_height,
                    units="in",
                    device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     diamide_path = snakemake@input[["diamide"]],
     aminoacid_path = snakemake@input[["aminoacid"]],
     nitrogen_path = snakemake@input[["nitrogen"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     svg_out = snakemake@output[["svg"]],
     pdf_out = snakemake@output[["pdf"]],
     png_out = snakemake@output[["png"]],
     grob_out = snakemake@output[["grob"]])

