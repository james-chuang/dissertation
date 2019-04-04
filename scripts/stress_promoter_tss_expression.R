main = function(theme_spec,
                fonts_path,
                genic_path,
                intra_path,
                counts_path,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)
    library(ggrepel)

    ttf_import(fonts_path)
    loadfonts()

    df = read_tsv(genic_path) %>%
        select(tss_name,
               feature_name) %>%
        mutate(category="genic") %>%
        bind_rows(read_tsv(intra_path) %>%
                      select(tss_name,
                             feature_name) %>%
                      mutate(category="intragenic")) %>%
        left_join(read_tsv(counts_path) %>%
                      select(name, 9:ncol(.)) %>%
                      gather(key=sample, value=rlog_counts, -name),
                  by=c("tss_name"="name")) %>%
        group_by(tss_name,
                 feature_name,
                 category) %>%
        summarise(rlog_counts_mean = mean(rlog_counts),
                  rlog_counts_sd = sd(rlog_counts)) %>%
        group_by(category) %>%
        mutate(rank = rank(rlog_counts_mean),
               y = rank/max(rank))

    figure_2b = ggplot() +
        geom_segment(data = df,
                     aes(x=rlog_counts_mean-rlog_counts_sd,
                         xend=rlog_counts_mean+rlog_counts_sd,
                         y=y,
                         yend=y,
                         color=category),
                        alpha=0.4,
                        size=0.3) +
        geom_point(data = df,
                   aes(x=rlog_counts_mean,
                       y=y,
                       color=category),
                   shape=16,
                   alpha=0.95,
                   size=0.4) +
        geom_label(data = df %>%
                       filter(category=="intragenic"),
                   aes(x=rlog_counts_mean-rlog_counts_sd-0.1,
                       y=y,
                       label=feature_name),
                  hjust=1,
                  size=4/72*25.4,
                  label.r = unit(0, "pt"),
                  label.padding = unit(1, "pt"),
                  label.size = NA,
                  family="FreeSans") +
        scale_y_continuous(name="percentile",
                           expand = c(0,0.03),
                           breaks = scales::pretty_breaks(n=3),
                           labels = function(x) 100*x) +
        scale_x_continuous(name=expression(list(log[2] (normalized ~ counts), regularized)),
                           breaks = scales::pretty_breaks(n=4)) +
        ggtitle("TSS expression levels in oxidative stress",
                subtitle = "oxidative stress-induced promoters") +
        scale_color_tableau() +
        theme_default +
        theme(legend.position = c(0.8, 0.40),
              plot.subtitle = element_text(size=7))

    ggplot2::ggsave(pdf_out, plot=figure_2b, width=fig_width, height=fig_height, units="in")
    embed_fonts(pdf_out)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path = snakemake@input[["fonts_path"]],
     genic_path = snakemake@input[["genic_path"]],
     intra_path = snakemake@input[["intra_path"]],
     counts_path = snakemake@input[["counts_path"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

