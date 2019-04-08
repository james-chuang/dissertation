main = function(theme_spec,
                fonts_path,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)

    ttf_import(fonts_path)
    loadfonts()

    orf_label = textGrob(label="coding DNA sequence",
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

    ggplot2::ggsave(pdf_out, plot=diagram, width=fig_width, height=fig_height, units="in")
    embed_fonts(pdf_out)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path = snakemake@input[["fonts_path"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

