
main = function(theme_spec = "thesis_theme.R",
                data_path = "2018_11_07_cell_counts.tsv",
                fig_width,
                fig_height,
                frame_1,
                frame_2,
                frame_3){
    source(theme_spec)
    library(gganimate)

    df = read_tsv(data_path) %>%
        group_by(treatment_time, diamide, experiment_type) %>%
        summarise(ratio_mean=mean(ratio),
                  ratio_sd=sd(ratio),
                  proportion_mean=mean(proportion),
                  proportion_sd=sd(proportion)) %>%
        filter(diamide %in% c(0,0.9,1,1.25) &
                   experiment_type == "mutant vs. WT") %>%
        mutate(frame=case_when(diamide==0 ~ 1,
                               diamide %in% c(0.9, 1) ~ 2,
                               TRUE ~ 3))

    figure = ggplot(data=df) +
        geom_hline(yintercept=0.5,
                   linetype="dashed",
                   size = 1,
                   color="grey70") +
        geom_ribbon(aes(x=treatment_time,
                        ymax=proportion_mean+proportion_sd,
                        ymin=proportion_mean-proportion_sd,
                        group=diamide),
                    alpha=0.6) +
        geom_line(aes(x=treatment_time,
                      y=proportion_mean,
                      group=diamide,
                      color=diamide),
                  size = 1) +
        scale_x_continuous(expand=c(0,0),
                           name="days in diamide",
                           breaks=0:2) +
        scale_y_continuous(name=expression(atop("percent", italic("dsk2-pace") * " cells")),
                           breaks=c(0,0.25,0.5),
                           limits=c(0, max(df[["proportion_mean"]] +
                                               df[["proportion_sd"]]) * 1.05),
                           expand=c(0,0),
                           labels = function(x)100*x) +
        scale_color_viridis(breaks = c(0,0.9,1,1.25),
                            guide=guide_colorbar(barwidth=0.8,
                                                 barheight=7),
                            name="diamide (mM)",
                            option="viridis") +
        theme_default_presentation +
        theme(axis.text = element_text(size=12),
              axis.title.y = element_text(angle=0,
                                          vjust=0.5,
                                          margin=margin(r=3, unit="pt")),
              legend.position = "right",
              legend.justification = c(0.5, 0.5),
              legend.title = element_text(size=12,
                                          color="black",
                                          hjust=0),
              legend.text = element_text(size=10),
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank()) +
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
     data_path = snakemake@input[["data_path"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     frame_1 = snakemake@output[["frame_1"]],
     frame_2 = snakemake@output[["frame_2"]],
     frame_3 = snakemake@output[["frame_3"]])

