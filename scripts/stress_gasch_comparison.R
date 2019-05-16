import_gasch = function(path){
    read_tsv(path) %>%
        gather(condition, lfc, -name) %>%
        mutate_at(vars(name, condition), ~fct_inorder(., ordered=TRUE)) %>%
        separate(condition, into=c("stress", "timepoint")) %>%
        return()
}

join_with_gasch = function(path, gasch, stress_id, timepoint_id){
    read_tsv(path) %>%
        group_by(genic_name) %>%
        summarise(tfiib_lfc = mean(log2_foldchange, na.rm=TRUE)) %>%
        inner_join(gasch %>%
                       filter(stress==stress_id & timepoint==timepoint_id) %>%
                       select(-c(stress, timepoint), gasch_lfc=lfc),
                   by=c("genic_name"="name"))
}

plot_scatter = function(df, xtitle, ytitle){
    plot = ggplot(data = df,
           aes(x=gasch_lfc, y=tfiib_lfc)) +
        geom_vline(xintercept = 0,
                   size=0.2,
                   color="grey70") +
        geom_hline(yintercept = 0,
                   size=0.2,
                   color="grey70") +
        geom_abline(slope=1, intercept = 0,
                    size=0.2,
                    color="grey70") +
        stat_bin_hex(geom="point",
                     aes(color=..count..),
                     binwidth=c(0.13, 0.13),
                     size=0.3,
                     shape=16,
                     alpha=0.8) +
        annotate(geom="label",
                 x=min(df[["gasch_lfc"]], na.rm=TRUE),
                 y=max(df[["tfiib_lfc"]], na.rm=TRUE),
                 label=paste0("R=",  cor(df[["tfiib_lfc"]],
                           df[["gasch_lfc"]],
                           use="complete.obs") %>%
                     signif(2)),
                 hjust=0,
                 vjust=1,
                 size=8/72*25.4,
                 label.size=NA,
                 label.padding=unit(1, "pt"),
                 label.r=unit(0, "pt"),
                 family="FreeSans") +
        scale_x_continuous(breaks=scales::pretty_breaks(n=4),
                           name=xtitle) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=4),
                           name=ytitle) +
        scale_color_viridis() +
        theme_default +
        theme(legend.position="none")
        # theme(legend.position="none",
        #       axis.title.y=element_text(angle=0,
        #                                 vjust=0.5))
    return(plot)
}

main = function(theme_spec = "thesis_theme.R",
                fonts_path,
                gasch_one_path = "gasch00_fig1data_cleaned.tsv",
                gasch_three_path = "gasch00_fig3data_cleaned.tsv",
                tfiib_diamide_path = "diamide-v-YPD_tfiib-chipnexus-libsizenorm-peaks-diffbind-results-genic-all.tsv",
                tfiib_aminoacid_path = "SD-v-SC_tfiib-chipnexus-libsizenorm-peaks-diffbind-results-genic-all.tsv",
                nitrogen_ontology_path,
                fig_width,
                fig_height,
                pdf_out){
    source(theme_spec)
    library(cowplot)

    ttf_import(fonts_path)
    loadfonts()

    gasch_one = import_gasch(gasch_one_path)
    gasch_three = import_gasch(gasch_three_path)

    diamide = join_with_gasch(path = tfiib_diamide_path,
                    gasch = gasch_one,
                    stress_id = "diamide",
                    timepoint_id = "40m")
    aminoacid = join_with_gasch(path = tfiib_aminoacid_path,
                                gasch = gasch_one,
                                stress_id = "aminoacid",
                                timepoint_id = "30m")

    scatter_diamide = plot_scatter(df = diamide,
                                   # xtitle = expression("RNA microarray log"[2] ~ textstyle(frac("45 min. oxidative stress", "YPD"))),
                                   xtitle = expression("RNA microarray log"[2] ~ frac("45 min. oxidative stress", "YPD")),
                                   ytitle = expression(atop("TFIIB ChIP-nexus",
                                                            "log"[2] ~ frac("40 min. oxidative stress", "YPD")))) +
        theme(plot.margin=margin(b=8, r=0, unit="pt"))
    scatter_aminoacid = plot_scatter(df = aminoacid,
                                     xtitle = expression("RNA microarray log"[2] ~ frac("30 min. amino acid stress", "SC")),
                                     ytitle = expression(atop("TFIIB ChIP-nexus",
                                                              "log"[2] ~ frac("30 min. amino acid stress", "SC")))) +
        theme(plot.margin=margin(t=8, r=0, unit="pt"))

    figure = plot_grid(scatter_diamide,
                       scatter_aminoacid,
                       ncol=1,
                       align="v",
                       axis="lr")

    ggplot2::ggsave(pdf_out,
                    plot=figure,
                    width=fig_width,
                    height=fig_height,
                    units="in",
                    device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path = snakemake@input[["fonts_path"]],
     gasch_one_path = snakemake@input[["gasch_one"]],
     gasch_three_path = snakemake@input[["gasch_three"]],
     tfiib_diamide_path = snakemake@input[["tfiib_diamide"]],
     tfiib_aminoacid_path = snakemake@input[["tfiib_aminoacid"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

