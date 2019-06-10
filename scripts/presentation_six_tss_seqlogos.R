
main = function(theme_spec,
                data_paths,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)
    library(ggseqlogo)
    #we don't use ggseqlogo plotting because it doesn't allow prior to be taken into account
    #we just use their dataframes with font information for geom_polygon

    slop=20
    tss_classes = c('genic', 'intragenic', 'antisense')

    df = tibble()
    for (i in 1:length(data_paths)){
        df = read_tsv(data_paths[i],
                      comment="#",
                      col_names=c('position','A','C','G','T','entropy','low','high','weight')) %>%
            mutate(position = position - as.integer(slop+1)) %>%
            gather(key=base, value=count, c('A','C','G','T')) %>%
            group_by(position) %>%
            mutate(height=entropy*count/sum(count)) %>%
            arrange(height, .by_group=TRUE) %>%
            mutate(base_low = lag(cumsum(height), default=0),
                   base_high = base_low+height) %>%
            left_join(ggseqlogo:::get_font("helvetica_bold"), by=c("base"="letter")) %>%
            group_by(position, base) %>%
            mutate(x = scales::rescale(x, to=c((1-(first(weight)))/2, 1-(1-first(weight))/2)),
                   x = x+(position-0.5),
                   y = scales::rescale(y, to=c(first(base_low), first(base_high))),
                   tss_class = tss_classes[i]) %>%
            bind_rows(df, .)
    }
    df %<>%
        mutate(tss_class = fct_inorder(tss_class, ordered=TRUE))

    fig_five_c = ggplot() +
        geom_polygon(data = df, aes(x=x, y=y, group=interaction(position, base), fill=base),
                     alpha=0.95) +
        geom_label(data=df %>%
                    distinct(tss_class),
                  aes(label=tss_class),
                  x=-12,
                  y=max(df[["y"]])*0.8,
                  size=12/72*25.4,
                  label.size=NA,
                  label.padding=unit(2, "pt"),
                  label.r = unit(0, "pt"),
                  hjust=0,
                  family="FreeSans") +
        scale_fill_manual(values = c('#109648', '#255C99', '#F7B32B', '#D62839', '#D62839'),
                          breaks = c('A','C','G','T','U')) +
        scale_y_continuous(limits = c(NA, max(df[["y"]]) *1.05 ),
                           breaks = c(0, 0.4),
                           expand=c(0,0),
                           name = "bits") +
        scale_x_continuous(limits = c(-12.5, 12.5),
                           expand = c(0,0),
                           labels = function(x)case_when(x==0 ~ "TSS",
                                                         x==10 ~ "+10 nt",
                                                         x>0 ~ paste0("+", x),
                                                         TRUE ~ as.character(x))) +
        facet_grid(tss_class~., switch="y") +
        ggtitle(bquote(italic("spt6-1004") ~ "TSSs")) +
        theme_default_presentation +
        theme(legend.position = "none",
              axis.title.y = element_text(angle=0, hjust=1, vjust=0.5),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size=12, color="black", face="plain",
                                         margin=margin(1,0,0,0,"pt")),
              axis.line = element_line(size=0.25, color="grey65"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              plot.margin = margin(1,0,0,0,"pt"))

    ggsave(pdf_out,
           plot=fig_five_c,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     data_paths = snakemake@input[["data_paths"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

