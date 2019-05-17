
main = function(theme_spec = "thesis_theme.R",
                data_path = "coding-transcripts-nonoverlapping-scaled_RNA-seq-sense.tsv.gz",
                fig_width = 3,
                fig_height = 2,
                pdf_out = "test.pdf"){
    source(theme_spec)

    df = read_tsv(data_path,
                  col_names = c("group",
                                "sample",
                                "annotation",
                                "assay",
                                "index",
                                "position",
                                "signal")) %>%
        group_by(group,
                 index,
                 position) %>%
        summarise(signal = mean(signal, na.rm=TRUE)) %>%
        group_by(index) %>%
        mutate(signal = (signal-mean(signal, na.rm=TRUE))/sd(signal, na.rm=TRUE)) %>%
        group_by(group,
                 position) %>%
        summarise(low = quantile(signal, 0.25, na.rm=TRUE),
                  mid = median(signal, na.rm=TRUE),
                  high = quantile(signal, 0.75, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(group = ordered(group,
                               levels=c("non-depleted",
                                        "depleted")))

    figure = ggplot(data = df,
           aes(x=position,
               ymin=low,
               ymax=high,
               y=mid,
               fill=group,
               color=group)) +
        geom_vline(xintercept = c(0,2),
                   color = "grey70",
                   size = 0.3) +
        geom_ribbon(alpha=0.2,
                    linetype="blank") +
        geom_line(alpha=0.85) +
        scale_x_continuous(expand = c(0,0),
                           breaks = c(0,2),
                           labels = c("TSS",
                                      "CPS"),
                           name = NULL) +
        scale_y_continuous(breaks = scales::pretty_breaks(n=2),
                           name = "standard score",
                           limits = c(min(df[["low"]]),
                                      max(df[["high"]]) * 1.05),
                           expand = c(0,0)) +
        scale_fill_few() +
        scale_color_few() +
        ggtitle("sense RNA-seq signal") +
        theme_default +
        theme(legend.position = c(0.55, 0.85),
              legend.key.height = unit(12, "pt"),
              axis.text.x = element_text(size=10),
              panel.grid = element_blank(),
              plot.margin = margin(t=1, unit="pt"))

    ggsave(pdf_out,
           plot=figure,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     data_path = snakemake@input[["data_path"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

