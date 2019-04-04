import = function(path, fdr_cutoff_tss, category_id){
    read_tsv(path) %>%
        filter(tfiib_start >= 0 & ! is.na(tss_FDR)) %>%
        select(tss_name, tss_lfc, tss_FDR, tfiib_name, tfiib_lfc, tfiib_FDR) %>%
        distinct(tss_name, .keep_all = TRUE) %>%
        mutate(significance = (tss_FDR > -log10(fdr_cutoff_tss)),
               direction = sign(tss_lfc)) %>%
        count(significance, direction) %>%
        mutate(category = category_id)
}

main = function(theme_spec,
                fonts_path,
                fdr_cutoff_tss,
                genic_path,
                intra_path,
                anti_path,
                inter_path,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)

    ttf_import(fonts_path)
    loadfonts()

    df = import(path=genic_path, fdr_cutoff_tss=fdr_cutoff_tss, category_id="genic") %>%
        bind_rows(import(path=intra_path, fdr_cutoff_tss=fdr_cutoff_tss, category_id="intragenic")) %>%
        bind_rows(import(path=anti_path, fdr_cutoff_tss=fdr_cutoff_tss, category_id="antisense")) %>%
        bind_rows(import(path=inter_path, fdr_cutoff_tss=fdr_cutoff_tss, category_id="intergenic")) %>%
        filter(significance) %>%
        mutate(category = fct_inorder(category, ordered=TRUE) %>%
                   fct_rev())

    diffexp_summary = ggplot() +
        geom_col(data = df,
                 aes(x=category, y=n*direction, fill=factor(direction)),
                 color="white",
                 width=1,
                 alpha=0.7) +
        geom_segment(data = df %>%
                         filter(direction==-1 &
                                   category != "genic"),
                     aes(x=category, xend=category,
                         yend=-n),
                     y=-500,
                     position = position_nudge(x=-0.17),
                     alpha=0.2, size=0.2) +
        geom_text(data = df,
                  aes(x=category,
                      y=if_else(n>40,
                                n*direction/2,
                                n*direction+5*direction),
                      label=n,
                      hjust=case_when(direction==1 & n<40 ~ 0,
                                      direction==-1 & n<40 ~ 1,
                                      TRUE ~ 0.5)),
                  size=9/72*25.4,
                  family="FreeSans") +
        geom_text(data = df %>% distinct(category),
                  aes(x=category, y=-500, label=category),
                  hjust=0, vjust=0.67,
                  size=7/72*25.4,
                  family="FreeSans") +
        annotate(geom="label",
                 x=4.8, y=175/2, label="upregulated",
                 fill = "#af8dc3",
                 label.r = unit(0, "pt"),
                 label.size = NA,
                 size = 7/72*25.4,
                 family="FreeSans",
                 alpha=0.7) +
        annotate(geom="label",
                 x=4.8, y=-551/2, label="downregulated",
                 fill = "#7fbf7b",
                 label.r = unit(0, "pt"),
                 label.size = NA,
                 size = 7/72*25.4,
                 family="FreeSans",
                 alpha=0.7) +
        scale_x_discrete(expand = c(0,1)) +
        coord_flip() +
        scale_fill_manual(values = c("#7fbf7b", "#af8dc3"),
                          guide=FALSE) +
        theme_void() +
        theme(plot.margin = margin(b=-10, unit="pt"))

    orf_label = textGrob(label="ORF",
                         x=0.55, y=0.5,
                         gp=gpar(fontsize=7,
                                 fontfamily="FreeSans"))
    orf_box = polygonGrob(x=c(0.375, 0.725*0.93, 0.725, 0.725*0.93, 0.375),
                          y=c(0.6, 0.6, 0.5, 0.4, 0.4),
                          gp = gpar(fill="grey80",
                                    lwd=NA))
    inter = textGrob(label = "intergenic",
                     x=0.06, y=0.35, hjust=0, vjust=0.5,
                     gp=gpar(fontsize=7, fontfamily="FreeSans"))
    inter_dash = linesGrob(x=c(0.31, 0.31), y=c(0.45, 0.35), gp=gpar(lty="twodash"))
    inter_arrow = linesGrob(x=c(0.31, 0.245), y=c(0.35, 0.35),
                            arrow = arrow(length=unit(0.15, "cm")))
    intra = textGrob(label = "intragenic",
                     x=0.78, y=0.93, hjust=0, vjust=0.5,
                     gp=gpar(fontsize=7, fontfamily="FreeSans"))
    intra_dash = linesGrob(x=c(0.6, 0.6), y=c(0.65, 0.93), gp=gpar(lty="twodash"))
    intra_arrow = linesGrob(x=c(0.6, 0.77), y=c(0.93, 0.93),
                            arrow = arrow(length=unit(0.15, "cm")))
    genic = textGrob(label = "genic",
                     x=0.78, y=0.73, hjust=0, vjust=0.5,
                     gp=gpar(fontsize=7, fontfamily="FreeSans"))
    genic_dash = linesGrob(x=c(0.35, 0.35), y=c(0.56, 0.73), gp=gpar(lty="twodash"))
    genic_arrow = linesGrob(x=c(0.35, 0.77), y=c(0.73, 0.73),
                            arrow = arrow(length=unit(0.15, "cm")))
    anti = textGrob(label = "antisense",
                     x=0.19, y=0.2, hjust=0, vjust=0.5,
                     gp=gpar(fontsize=7, fontfamily="FreeSans"))
    anti_dash = linesGrob(x=c(0.5, 0.5), y=c(0.35, 0.2), gp=gpar(lty="twodash"))
    anti_arrow = linesGrob(x=c(0.5, 0.37), y=c(0.2, 0.2),
                            arrow = arrow(length=unit(0.15, "cm")))
    genome_line = linesGrob(x=c(0.09,0.93), y=c(0.5, 0.5),
                            gp=gpar(col="grey60"))
    transcript_line = linesGrob(x=c(0.35,0.77), y=c(0.5, 0.5),
                                gp=gpar(lwd=2,
                                        col="grey25",
                                        lineend="square"))

    diagram = gTree(children = gList(genome_line,
                                     transcript_line,
                                     intra_dash, intra_arrow,
                                     genic_dash, genic_arrow,
                                     inter_dash, inter_arrow,
                                     anti_dash, anti_arrow,
                                     orf_box, orf_label, inter, intra, anti, genic))

    figure_2a = arrangeGrob(diagram, diffexp_summary, ncol=1, heights=c(0.5, 1))

    ggplot2::ggsave(pdf_out, plot=figure_2a, width=fig_width, height=fig_height, units="in")
    embed_fonts(pdf_out)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path = snakemake@input[["fonts_path"]],
     fdr_cutoff_tss = snakemake@params[["fdr_cutoff_tss"]],
     genic_path = snakemake@input[["genic_path"]],
     intra_path = snakemake@input[["intra_path"]],
     anti_path = snakemake@input[["anti_path"]],
     inter_path = snakemake@input[["inter_path"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

