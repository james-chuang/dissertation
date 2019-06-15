
main = function(theme_spec = "thesis_theme.R",
                tfiib_data = "SWC4_all-assays.tsv.gz",
                fig_width = 14,
                fig_height = 6.65,
                frame_1 = "test/wasda_test0001.svg",
                frame_2 = "test/wasda_test0002.svg"){
    source(theme_spec)
    library(gganimate)

    sample_list = c("WT-37C-1", "WT-37C-2", "spt6-1004-37C-1", "spt6-1004-37C-2")
    transcript_end = 1.676
    orf_start = 0.009
    orf_end = 1.440
    gene_id = "SWC4"

    df = read_tsv(tfiib_data,
                  col_names = c("group", "sample", "annotation", "assay", "index", "position", "signal")) %>%
        filter(sample %in% sample_list &
                   assay=="TFIIB-ChIP-nexus-protection") %>%
        group_by(group, position) %>%
        summarise(signal = mean(signal)) %>%
        ungroup() %>%
        mutate(group = ordered(group,
                               levels = c("spt6+", "spt6-1004-37C"),
                               labels = c("wild-type", "spt6-1004"))) %>%
        filter(between(position, -0.35, 2.2))

    df %<>% mutate(frame = as.integer(group))

    max_signal = max(df[["signal"]])
    orf_height = max_signal*0.1
    anno_y = max_signal*0.93
    min_position = min(df[["position"]])
    max_position = max(df[["position"]])

    figure = ggplot(data=df,
                       aes(x=position,
                           y=signal,
                           fill=group,
                           group=group)) +
        geom_area(alpha=0.8,
                  position=position_identity()) +
        scale_x_continuous(expand=c(0,0),
                           breaks=scales::pretty_breaks(3),
                           labels=function(x) case_when(x == 0 ~ "TSS",
                                                        x %>% near(2) ~ paste(x, "kb"),
                                                        TRUE ~ as.character(x))) +
        scale_y_continuous(oob=scales::squish,
                           expand=c(0,0),
                           breaks=scales::pretty_breaks(n=2),
                           name="normalized counts") +
        annotate(geom="segment",
                 x=0,
                 xend=transcript_end,
                 y=anno_y,
                 yend=anno_y,
                 color="grey50") +
        annotate(geom="segment",
                 x=c(min_position, 1.863),
                 xend=c(-0.118, max_position),
                 y=anno_y,
                 yend=anno_y,
                 color="grey50") +
        annotate(geom="polygon",
                 x=c(orf_start,
                     (orf_end-orf_start)*.93+orf_start,
                     orf_end,
                     (orf_end-orf_start)*.93+orf_start,
                     orf_start),
                 y=c(anno_y+orf_height/2,
                     anno_y+orf_height/2,
                     anno_y,
                     anno_y-orf_height/2,
                     anno_y-orf_height/2),
                 fill="grey80") +
        annotate(geom="rect",
                 xmin=c(min_position, 1.881),
                 xmax=c(-0.216, max_position),
                 ymin=anno_y-orf_height/2,
                 ymax=anno_y+orf_height/2,
                 fill="grey80") +
        annotate(geom="text",
                 label=gene_id,
                 x=(orf_end-orf_start)/2+orf_start,
                 y=anno_y,
                 vjust = 0.5,
                 size=12/72*25.4,
                 family="FreeSans",
                 fontface="italic") +
        scale_fill_ptol(labels=c("wild-type",
                                 bquote(italic("spt6-1004")))) +
        ggtitle("TFIIB protection") +
        theme_default_presentation +
        theme(legend.key.height = unit(14, "pt"),
              legend.key.width = unit(25, "pt"),
              legend.spacing.y = unit(0, "pt"),
              legend.position = c(0.7, 0.7),
              axis.title.y = element_text(margin=margin(r=1, unit="pt")),
              axis.text.x = element_text(size=12),
              axis.title.x = element_blank(),
              panel.border = element_blank(),
              axis.line.y = element_line(size=0.25, color="grey70"),
              panel.grid = element_blank()) +
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
     tfiib_data = snakemake@input[["tfiib_data"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     frame_1 = snakemake@output[["frame_1"]],
     frame_2 = snakemake@output[["frame_2"]])

