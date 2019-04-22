
import = function(path,
                  sample_list,
                  assay_name){
    read_tsv(path, col_names = c("group", "sample", "annotation", "assay", "index", "position", "signal")) %>%
        filter((sample %in% sample_list) & ! is.na(signal)) %>%
        mutate(group = fct_inorder(group, ordered=TRUE)) %>%
        group_by(group, index, position) %>%
        summarise(signal=mean(signal)) %>%
        group_by(index) %>%
        mutate(signal=(signal-mean(signal))/sd(signal)) %>%
        group_by(group, position) %>%
        summarise(mid=median(signal, na.rm=TRUE),
                  low=quantile(signal, 0.25, na.rm=TRUE),
                  high=quantile(signal, 0.75, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(assay=assay_name) %>%
        return()
}

main = function(theme_spec,
                fonts_path,
                spt5_data,
                netseq_data,
                rnapii_data,
                annotation,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)

    ttf_import(fonts_path)
    loadfonts()

    sample_list = c("non-depleted-1", "non-depleted-2", "depleted-1", "depleted-2")

    df = import(netseq_data,
                sample_list=c("non-depleted-3",
                              "non-depleted-4",
                              "depleted-1",
                              "depleted-2"),
                assay_name="sense NET-seq") %>%
        bind_rows(import(rnapii_data,
                    sample_list=c("non-depleted-1",
                                  "non-depleted-2",
                                  "depleted-1",
                                  "depleted-2"),
                    assay_name="RNAPII ChIP-seq")) %>%
        bind_rows(import(spt5_data,
                    sample_list=c("non-depleted-1",
                                  "non-depleted-2",
                                  "depleted-1",
                                  "depleted-2"),
                    assay_name="Spt5 ChIP-seq")) %>%
        mutate(assay=fct_inorder(assay, ordered=TRUE) %>%
                        fct_rev())

    plot = ggplot(data = df, aes(x=position, color=group, fill=group)) +
        geom_vline(xintercept = c(0, 2), size=0.3, color="grey70") +
        geom_ribbon(aes(ymin=low, ymax=high), alpha=0.17, linetype='blank') +
        geom_line(aes(y=mid), alpha=0.75) +
        scale_x_continuous(breaks = c(0, 2),
                           labels = c("TSS", "CPS"),
                           name = NULL,
                           expand = c(0,0)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n=3),
                           name = "standard score") +
        facet_wrap(~assay, ncol=1, scales="free_y") +
        scale_color_ptol() +
        scale_fill_ptol() +
        theme_default +
        theme(legend.key.height = unit(10, "pt"),
              strip.text=element_text(size=9, hjust=0),
              panel.grid = element_blank(),
              legend.position = c(0.7, 0.27),
              legend.background = element_rect(color=NA, fill="white", size=0))

    ggsave(pdf_out,
           plot=plot,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path=snakemake@input[["fonts_path"]],
     spt5_data = snakemake@input[["spt5_data"]],
     netseq_data = snakemake@input[["netseq_data"]],
     rnapii_data = snakemake@input[["rnapii_data"]],
     annotation = snakemake@input[["annotation"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

