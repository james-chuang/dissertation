import = function(path){
    read_tsv(path,
             col_names=c("group",
                         "factor",
                         "annotation",
                         "index",
                         "position",
                         "signal")) %>%
        return()
}

main = function(theme_spec,
                data_paths,
                annotation,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)
    library(ggridges)

    df=tibble()
    for (path in data_paths){
        df %<>%
            bind_rows(import(path))
    }

    df %<>%
        group_by(factor) %>%
        mutate(signal=scales::rescale(signal)) %>%
        spread(factor, signal) %>%
        mutate_at(vars(Ser2P, Ser5P),
                  ~log2(./Rpb1)) %>%
        group_by(index) %>%
        mutate(Rpb1=(Rpb1 - mean(Rpb1, na.rm=TRUE)) / sd(Rpb1, na.rm=TRUE)) %>%
        gather(factor, enrichment, c(Ser2P, Ser5P)) %>%
        group_by(factor) %>%
        mutate(enrichment=pmin(enrichment, max(enrichment[!is.infinite(enrichment)])),
               enrichment=pmax(enrichment, min(enrichment[!is.infinite(enrichment)])),
               enrichment=(enrichment-mean(enrichment,na.rm=TRUE))/sd(enrichment, na.rm=TRUE)) %>%
        group_by(group, position, factor) %>%
        summarise(rnap=median(Rpb1, na.rm=TRUE),
                  enrichment=median(enrichment, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(rnap=scales::rescale(rnap),
               group=ordered(group, levels=c("non-depleted",
                                             "depleted")))

    enrichment_plot = ggplot() +
        geom_vline(xintercept=0,
                   color="grey70",
                   size=0.5) +
        geom_density_ridges_gradient(data=df,
                                aes(x=position,
                                    y=factor,
                                    height=rnap,
                                    fill=enrichment),
                                stat="identity",
                                size=0.3,
                                scale=0.9,
                                panel_scaling=FALSE) +
        geom_hline(yintercept=c(1,2),
                   size=0.3) +
        geom_label(data = df %>%
                       filter(group=="non-depleted") %>%
                       distinct(group, factor),
                   aes(label=paste0(factor, "\nenrichment"),
                       y=as.numeric(as.factor(factor))+1),
                   x=-0.4,
                   hjust=0,
                   vjust=1.1,
                   size=10/72*25.4,
                   label.r=unit(0, "pt"),
                   label.size=NA,
                   label.padding=unit(2,"pt"),
                   alpha=0.6,
                   family="FreeSans") +
        facet_grid(.~group) +
        scale_fill_viridis(guide=guide_colorbar(title="relative\nenrichment",
                                                title.hjust=0,
                                                barheight=4.9,
                                                barwidth=0.5),
                           option="viridis") +
        scale_x_continuous(expand=c(0,0),
                           breaks=scales::pretty_breaks(n=3),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        x==1.5 ~ paste(x, "kb"),
                                                        TRUE ~ as.character(x))) +
        scale_y_discrete(expand=c(0,0),
                         name="RNA Pol II\nlevels") +
        theme_default +
        theme(axis.title.x=element_blank(),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5,
                                        hjust=1),
              axis.text.x=element_text(size=10),
              axis.text.y=element_blank(),
              strip.text=element_text(color="black",
                                      size=10,
                                      hjust=0),
              legend.position="right",
              legend.justification=c(0, 0.5),
              legend.title=element_text(size=10),
              legend.text=element_text(size=8),
              panel.grid = element_blank())

    ggsave(pdf_out,
           plot=enrichment_plot,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     data_paths = snakemake@input[["data_paths"]],
     annotation = snakemake@input[["annotation"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

