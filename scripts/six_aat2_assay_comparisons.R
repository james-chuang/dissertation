
main = function(theme_spec,
                data_path,
                fonts_path,
                pdf_out, fig_height, fig_width){
    source(theme_spec)
    library(cowplot)

    ttf_import(fonts_path)
    loadfonts()

    df = read_tsv(data_path,
                  col_names = c('condition', 'sample', 'annotation',
                                'assay', 'index', 'position', 'signal')) %>%
        filter(condition == "spt6-1004-37C") %>%
        group_by(assay, position) %>%
        summarise(signal = mean(signal)) %>%
        ungroup() %>%
        mutate(assay = ordered(assay,
                               levels = c("RNA-seq-sense",
                                          "TSS-seq-sense",
                                          "TFIIB-ChIP-nexus-protection"),
                               labels = c("\"RNA-seq (Uwimana \" * italic(\"et al.\") * \", 2017)\"",
                                          "\"TSS-seq\"",
                                          "\"TFIIB ChIP-nexus protection\"")))

    coverage = ggplot(data = df %>%
               filter(position %>% between(-0.2, 1.6)),
           aes(x=position, y=signal)) +
        geom_col(size=0.3, color=viridis(1)) +
        facet_wrap(~assay,
                   scales= "free_y",
                   ncol=1,
                   labeller = label_parsed) +
        scale_x_continuous(expand=c(0,0),
                           labels = function(x)case_when(x==0 ~ "TSS",
                                                         x==1.5 ~ paste(x, "kb"),
                                                         TRUE ~ as.character(x))) +
        scale_y_continuous(name = "normalized signal",
                           limits = c(0, NA),
                           expand = c(0,0),
                           breaks = scales::pretty_breaks(n=2)) +
        theme_default +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size=10, family="FreeSans"),
              strip.text = element_text(hjust=0, family="FreeSans"),
              panel.spacing.y = unit(4, "pt"),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              axis.line.y = element_line(size=0.2,
                                         color="black"),
              axis.ticks = element_line(color="black"),
              plot.margin = margin(t=-0.3, unit="cm"))

    diagram = ggplot() +
        annotate(geom="segment", color="grey50",
                 x=0, xend=1.348, y=0, yend=0) +
        annotate(geom="polygon", fill="grey80",
                 x=c(.028, (1.257)*.93+0.028, 1.257,
                     (1.257)*.93+0.028, .028),
                 y=c(1,1,0,-1,-1)) +
        annotate(geom="text",
                 label="AAT2",
                 family="FreeSans",
                 fontface="italic",
                 x=(1.257)/2+0.028, y=0,
                 size=10/72*25.4) +
        scale_x_continuous(limits = c(-0.2, 1.6),
                           expand=c(0,0)) +
        scale_y_continuous(limits = c(-1, 1), expand=c(0,0)) +
        theme_void() +
        theme(plot.margin = margin(-2,0,-0.3,0,"cm"))

    plot = plot_grid(diagram, coverage,
                     ncol = 1,
                     rel_heights = c(0.17,1),
                     rel_widths = c(1,1),
                     align = "vh",
                     axis = "trbl")

    ggplot2::ggsave(pdf_out,
                    plot=plot,
                    height=fig_height,
                    width=fig_width,
                    units="in",
                    device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     data_path = snakemake@input[["data_path"]],
     fonts_path = snakemake@input[["fonts_path"]],
     pdf_out = snakemake@output[["pdf"]],
     fig_height = snakemake@params[["height"]],
     fig_width = snakemake@params[["width"]])

