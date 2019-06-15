
import = function(path,
                  normalize=FALSE){
    df = read_tsv(path,
                  col_names=c('group', 'sample', 'annotation',
                              'assay', 'index', 'position', 'signal')) %>%
        group_by(annotation) %>%
        mutate(anno_group = if_else(str_detect(annotation, "intragenic"),
                                    "intragenic",
                                    "wild-type genic"))
    if (normalize){
        df %<>%
            group_by(group, sample, annotation,
                     index) %>%
            mutate(signal = scales::rescale(signal))
    }
    df %>%
        group_by(group, annotation, anno_group, position) %>%
        summarise(mid = median(signal),
                  low = quantile(signal, 0.25),
                  high = quantile(signal, 0.75),
                  mean = mean(signal)) %>%
        ungroup() %>%
        mutate(group = ordered(group,
                               levels=c("spt6+", "spt6-1004-37C"),
                               labels=c("WT", "spt6-1004-37C")),
               anno_group = ordered(anno_group,
                                    levels=c("wild-type genic",
                                             "intragenic")),
               index = as.integer(fct_inorder(annotation))) %>%
        return()
}

main = function(theme_spec = "thesis_theme.R",
                mnase_data = "spt6-induced-intragenic-TSSs-MNase-clusters-v-genic_MNase-seq.tsv.gz",
                gc_data = "spt6-induced-intragenic-TSSs-MNase-clusters-v-genic_GC-pct.tsv.gz",
                fig_width=14,
                fig_height=5.65,
                pdf_out){
    source(theme_spec)
    library(cowplot)

    gc_df = import(gc_data) %>%
        filter(group=="WT" & annotation != "all intragenic")

    rr = max(gc_df[["mean"]]) - min(gc_df[["mean"]])

    gc_plot = ggplot(data=gc_df) +
        geom_vline(xintercept=0,
                   size=0.4,
                   color="grey70") +
        annotate(geom="rect",
                 xmin=-0.09,
                 xmax=-0.04,
                 ymin=min(gc_df[["mean"]]) - 0.05 * rr,
                 ymax=max(gc_df[["mean"]]) + 0.05 * rr,
                 fill="grey95") +
        geom_line(aes(x=position,
                      y=mean,
                      group=annotation,
                      color=anno_group),
                  size=0.75) +
        scale_x_continuous(expand=c(0,0),
                           breaks=scales::pretty_breaks(n=5),
                           labels = function(x){case_when(x==0 ~ "TSS",
                                                          near(x, 0.4) ~ "0.4 kb",
                                                          TRUE ~ as.character(x))},
                           name=NULL) +
        scale_y_continuous(expand=c(0,0),
                           name="mean GC%") +
        scale_color_viridis_d(option="viridis",
                              end=0.85) +
        theme_default_presentation +
        theme(panel.grid=element_blank(),
              axis.text.x = element_text(size=12),
              axis.title.y = element_text(angle=0,
                                          vjust=0.5,
                                          hjust=1),
              legend.justification=c(0,0),
              legend.position=c(0.02,0.2))

    ggplot2::ggsave(pdf_out,
                    plot=gc_plot,
                    width=fig_width,
                    height=fig_height,
                    units="cm",
                    device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     gc_data = snakemake@input[["gc_data"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

