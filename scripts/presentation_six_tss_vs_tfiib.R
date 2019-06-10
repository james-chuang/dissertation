
import = function(df, path, category){
    df = read_tsv(path) %>%
        mutate(category=category) %>%
        bind_rows(df, .) %>%
        return()
}

main = function(theme_spec,
                genic,
                intragenic,
                antisense,
                fig_width,
                fig_height,
                pdf_out){
    source(theme_spec)

    df = tibble() %>%
        import(genic, 'genic') %>%
        import(intragenic, 'intragenic') %>%
        import(antisense, 'antisense') %>%
        mutate(category=fct_inorder(category, ordered=TRUE)) %>%
        filter(! is.na(tss_fdr))
    count_df = df %>%
        count(category)

    fig_two_d = ggplot(data = df) +
        geom_hline(yintercept=0,
                   color="grey70",
                   size=0.5) +
        geom_vline(xintercept=0,
                   color="grey70",
                   size=0.5) +
        geom_abline(intercept=0,
                    slope=1,
                    color="grey70",
                    size=0.5) +
        stat_bin_hex(geom="point",
                     aes(x=tss_lfc, y=tfiib_lfc, color=(..count..)),
                     binwidth=c(0.14, 0.14),
                     alpha=0.8,
                     size=0.4,
                     shape=16,
                     fill=NA) +
        geom_label(data=count_df,
                  aes(label=paste0("n=",n)),
                  x=-6,
                  y=5,
                  hjust=0,
                  size=8/72*25.4,
                  label.padding = unit(2, "pt"),
                  label.r = unit(0, "pt"),
                  label.size = NA,
                  family="FreeSans") +
        facet_wrap(~category,
                   nrow=1) +
        scale_color_viridis(guide=FALSE,
                            option="inferno",
                            limits = c(NA, 20),
                            oob = scales::squish) +
        scale_y_continuous(limits = c(-4,NA),
                           name = expression(atop("TFIIB ChIP-nexus",
                                                  "log"[2] * textstyle(frac(italic("spt6-1004"), "wild-type"))))) +
        scale_x_continuous(name = expression("TSS-seq " * "log"[2] * textstyle(frac(italic("spt6-1004"), "wild-type")))) +
        theme_default_presentation +
        theme(strip.text = element_text(size=12, color="black",
                                        margin=margin(b=2, unit="pt")),
              axis.title.x = element_text(size=12),
              axis.title.y = element_text(size=12, angle=0, vjust=0.5),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major = element_line(size=0.1),
              panel.spacing.y = unit(0, "pt"))

    ggsave(pdf_out,
           plot=fig_two_d,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     genic = snakemake@input[["genic"]],
     intragenic = snakemake@input[["intragenic"]],
     antisense = snakemake@input[["antisense"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

