
import = function(path, sample_list){
    read_tsv(path, col_names = c("group", "sample", "annotation", "index", "position", "signal")) %>%
        filter((sample %in% sample_list) & ! is.na(signal)) %>%
        group_by(group, annotation, position) %>%
        summarise(mid = median(signal),
                  low = quantile(signal, 0.25),
                  high = quantile(signal, 0.75)) %>%
        ungroup() %>%
        mutate(group = ordered(group,
                               levels = c("WT-37C", "spt6-1004-37C"),
                               labels = c("WT", "italic(\"spt6-1004\")"))) %>%
        return()
}

main = function(theme_spec = "thesis_theme.R",
                mnase_data = "verified-transcripts-nonoverlapping-plusonenuc-allsamples-allannotations-midpoint-spikenorm.tsv.gz",
                fig_width = 14,
                fig_height = 6.65,
                frame_1 = "test/wasda_test0001.svg",
                frame_2 = "test/wasda_test0002.svg"){
    source(theme_spec)
    library(gganimate)

    sample_list = c("WT-37C-1", "spt6-1004-37C-1", "spt6-1004-37C-2")
    max_length=1.5

    df = import(mnase_data, sample_list=sample_list) %>%
        mutate(frame = as.integer(group))
        # mutate_at(vars(-c(group, annotation, position)), funs(.*10))

    figure = ggplot(data = df,
                       aes(x=position, color=group, fill=group)) +
        geom_ribbon(aes(ymin=low, ymax=high), alpha=0.3, linetype='blank') +
        geom_line(aes(y=mid),
                  size=1,
                  alpha=0.85) +
        scale_x_continuous(breaks = scales::pretty_breaks(n=3),
                           labels = function(x){case_when(x==0 ~ "+1 dyad",
                                                          x==max_length ~ paste(x, "kb"),
                                                          TRUE ~ as.character(x))},
                           name = NULL,
                           expand = c(0,0)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n=2),
                           limits = c(0, max(df[["high"]]*1.05)),
                           expand = c(0,0),
                           labels = function(x){if_else(x<0, abs(x), x)},
                           name = "normalized counts") +
        ggtitle("MNase-seq dyad signal",
                subtitle="3522 nonoverlapping coding genes") +
        scale_color_ptol(labels = c("wild-type", bquote(italic("spt6-1004")))) +
        scale_fill_ptol(labels = c("wild-type", bquote(italic("spt6-1004")))) +
        theme_default_presentation +
        theme(legend.key.height = unit(14, "pt"),
              legend.key.width = unit(25, "pt"),
              axis.text.x = element_text(size=12),
              legend.position = c(0.65, 0.80),
              panel.grid = element_blank(),
              plot.subtitle = element_text(margin=margin(0,0,-2.5,0,"pt")) ,
              plot.margin = margin(0, 15, 0, 0, "pt")) +
        transition_manual(frame,
                          cumulative=TRUE)

    animate(figure,
            renderer=file_renderer(dir = dirname(frame_1),
                                   str_extract(basename(frame_1), "[a-z,_]+"),
                                   overwrite=TRUE),
            device="svg",
            width=fig_width/2.54,
            height=fig_height/2.54)
}

main(theme_spec = snakemake@input[["theme"]],
     mnase_data = snakemake@input[["mnase_data"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     frame_1 = snakemake@output[["frame_1"]],
     frame_2 = snakemake@output[["frame_2"]])

