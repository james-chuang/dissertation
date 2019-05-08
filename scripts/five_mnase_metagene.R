
import = function(path, sample_list){
    read_tsv(path, col_names = c("group", "sample", "annotation", "index", "position", "signal")) %>%
        filter((sample %in% sample_list) & ! is.na(signal)) %>%
        mutate(group = fct_inorder(group, ordered=TRUE)) %>%
        group_by(group, annotation, position) %>%
        summarise(mid = median(signal),
                  low = quantile(signal, 0.25),
                  high = quantile(signal, 0.75)) %>%
        ungroup() %>%
        return()
}

main = function(theme_spec,
                fonts_path,
                mnase_data,
                fig_width,
                fig_height,
                pdf_out){
    source(theme_spec)

    ttf_import(fonts_path)
    loadfonts()

    sample_list = c("non-depleted-1, non-depleted-2", "non-depleted-3",
                    "depleted-1", "depleted-2", "depleted-3")

    df = import(mnase_data, sample_list=sample_list)

    plot = ggplot(data = df, aes(x=position, color=group, fill=group)) +
        geom_ribbon(aes(ymin=low, ymax=high), alpha=0.17, linetype='blank') +
        geom_line(aes(y=mid), alpha=0.75) +
        scale_x_continuous(breaks = seq(0, 1, 0.5),
                           labels = c("+1 dyad", "0.5", "1 kb"),
                           name = NULL,
                           expand = c(0,0)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n=2),
                           limits = c(0, max(df[["high"]])*1.05),
                           expand = c(0,0.002),
                           name = "normalized counts") +
        ggtitle("MNase-seq dyad signal") +
        scale_color_few() +
        scale_fill_few() +
        theme_default +
        theme(legend.key.height = unit(10, "pt"),
              panel.grid = element_blank(),
              legend.position = c(0.88, 0.95),
              legend.background = element_rect(color=NA, fill="white", size=0))

    ggsave(pdf_out,
           plot=plot,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path = snakemake@input[["fonts_path"]],
     mnase_data = snakemake@input[["mnase_data"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

