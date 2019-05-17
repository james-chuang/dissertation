
main = function(theme_spec,
                data_path,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)

    df = read_tsv(data_path) %>%
        mutate(group = ordered(group,
                               levels=c("WT", "spt6-1004")))

    som_plot = ggplot() +
        geom_vline(xintercept = 0,
                   color="grey70",
                   size=0.5) +
        geom_ribbon(data = df,
                    aes(ymin=low,
                        ymax=high,
                        x=position,
                        fill=group),
                    alpha=0.2) +
        geom_line(data = df,
                  aes(x=position,
                      y=mid,
                      color=group),
                  alpha=0.8) +
        geom_label(data = df %>%
                       distinct(unit, cluster, n),
                   aes(label=n,
                       fill=factor(cluster)),
                   x=max(df[["position"]])*0.92,
                   y=max(df[["high"]]),
                   hjust=1,
                   vjust=1,
                   label.size=NA,
                   label.r=unit(0,"pt"),
                   label.padding=unit(3,"pt"),
                   show.legend = FALSE,
                   size=8/72*25.4,
                   family="FreeSans") +
        facet_wrap(~unit,
                   ncol=8) +
        scale_x_continuous(expand=c(0,0),
                           breaks=scales::pretty_breaks(n=3),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        near(x,0.1) ~ paste(x, "kb"),
                                                        TRUE ~ as.character(x)),
                           name=NULL) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=2),
                           name="normalized dyad counts") +
        scale_color_ptol(labels=c("wild-type",
                                  expression(italic("spt6-1004")))) +
        scale_fill_manual(values=c("#7fc97f",
                                   "#beaed4",
                                   "#fdc086",
                                   rev(ptol_pal()(2))),
                          breaks=c(1,2,3),
                          labels=function(x) paste("cluster", x),
                          guide=guide_legend(override.aes=list(alpha=1))) +
        theme_default +
        theme(legend.position="right",
              legend.justification=c(1, 0.5),
              legend.spacing.x=unit(1,"pt"),
              legend.text.align=0,
              panel.grid=element_blank())

    ggsave(pdf_out,
           plot=som_plot,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     data_path = snakemake@input[["som_data"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

