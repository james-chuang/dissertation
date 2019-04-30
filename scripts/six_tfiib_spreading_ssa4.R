
main = function(theme_spec,
                fonts_path,
                tfiib_data,
                fig_width,
                fig_height,
                pdf_out){
    source(theme_spec)
    library(ggforce)

    ttf_import(fonts_path)
    loadfonts()

    sample_list = c("WT-37C-1", "WT-37C-2", "spt6-1004-37C-1", "spt6-1004-37C-2")

    df = read_tsv(tfiib_data,
                  col_names = c("group", "sample", "annotation", "assay", "index", "position", "signal")) %>%
        filter(sample %in% sample_list) %>%
        group_by(group, position) %>%
        summarise(signal = mean(signal)) %>%
        ungroup() %>%
        mutate(group = ordered(group,
                               levels = c("spt6+", "spt6-1004-37C"),
                               labels = c("WT", "spt6-1004")))

    fig_two_b = ggplot(data=df,
                       aes(x=position,
                           y=signal,
                           fill=group)) +
        geom_area(alpha=0.6,
                  position=position_identity()) +
        scale_x_continuous(expand=c(0,0),
                           breaks=scales::pretty_breaks(3),
                           labels=function(x) case_when(x == 0 ~ "TSS",
                                                        x %>% near(10) ~ "+10 kb",
                                                        x %>% near(2) ~ "+2 kb",
                                                        x > 0 ~ paste0("+", x),
                                                        TRUE ~ as.character(x))) +
        scale_y_continuous(limits = c(NA, 0.8),
                           oob=scales::squish,
                           expand=c(0,0),
                           breaks=c(0,0.8),
                           name="normalized counts") +
        facet_zoom(xlim=c(-0.3, 2.111+0.3),
                   zoom.size=1) +
        annotate(geom="segment",
                 x=0,
                 xend=2.111,
                 y=0.7,
                 yend=0.7,
                 color="grey50") +
        annotate(geom="polygon",
                 x=c(0.054,
                     1.929*.93+0.054,
                     1.983,
                     1.929*.93+0.054,
                     0.054),
                 y=c(0.75, 0.75, 0.7, 0.65, 0.65),
                 fill="grey80") +
        annotate(geom="text",
                 label="SSA4",
                 x=1.929/2+0.054,
                 y=0.7,
                 size=7/72*25.4,
                 family="FreeSans",
                 fontface="italic") +
        scale_fill_ptol(labels=c("wild-type",
                                 bquote(italic("spt6-1004")))) +
        ggtitle("TFIIB ChIP-nexus protection") +
        theme_default +
        theme(legend.key.height = unit(10, "pt"),
              legend.position = c(0.78, 0.85),
              strip.background = element_rect(fill="grey90", size=0, color=NA),
              axis.title.y = element_text(margin=margin(r=1, unit="pt")),
              axis.title.x = element_blank(),
              panel.border = element_blank(),
              axis.line.y = element_line(size=0.25, color="grey70"),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.spacing.y = unit(1, "pt"),
              plot.margin=margin(r=11, t=1, unit="pt"))

    ggplot2::ggsave(pdf_out,
                    plot=fig_two_b,
                    width=fig_width,
                    height=fig_height,
                    units="in",
                    device=cairo_pdf)
}


main(theme_spec = snakemake@input[["theme"]],
     fonts_path = snakemake@input[["fonts_path"]],
     tfiib_data = snakemake@input[["tfiib_data"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

