
main = function(theme_spec = "thesis_theme.R",
                blot_path = "spt6-western-3.png",
                fig_width = 7,
                fig_height = 3.25,
                pdf_out){
    source(theme_spec)
    library(png)

    text_edge = 0.17
    right_edge=0.10
    increment = (1-text_edge-right_edge)/4
    start = increment/2+0.02

    alignment_lines = segmentsGrob(x0=seq(0,3)*increment+start+text_edge,
                                   x1=seq(0,3)*increment+start+text_edge,
                                   y0=0,
                                   y1=1,
                                   gp = gpar(lty="dashed"))
    temp_lines = segmentsGrob(x0=c(0,2)*increment-0.06+text_edge+start,
                                 x1=c(1,3)*increment+0.06+text_edge+start,
                                 y0=0.86, y1=0.86)
    temp_labels = textGrob(label = c(paste0("30", "\U00B0", "C"),
                                     paste0("37", "\U00B0", "C")),
                              gp = gpar(fontsize=12,
                                        fontfamily="FreeSans"),
                              x=text_edge+start+c(1,5)*increment/2,
                              y=1,
                              vjust=1.1)
    antigen_labels = textGrob(label=c("Spt6",
                                      "Dst1",
                                      "Spt6:"),
                              x=text_edge,
                              y=c(0.53,
                                  0.17,
                                  0.78),
                              gp=gpar(fontsize=12,
                                      fontfamily="FreeSans"),
                              hjust=1)
    plus_labels = textGrob(label = "+",
                           x=text_edge+start+c(0,2)*increment,
                           y=0.78,
                           gp=gpar(fontsize=12,
                                   fontfamily="FreeSans"),
                           vjust=0.5)
    mutant_labels = textGrob(label = "1004",
                             x=text_edge+start+c(1,3)*increment,
                             y=0.78,
                             gp=gpar(fontsize=12,
                                     fontface="italic",
                                     fontfamily="FreeSans"),
                             vjust=0.5)
    blot = rasterGrob(readPNG(blot_path),
                      width=1-text_edge+0.04,
                      x=text_edge + (1-text_edge)/2,
                      height=0.85,
                      y=0.79,
                      vjust=1)
    spt6_outline = rectGrob(width=1-text_edge-0.06,
                            x=text_edge + (1-text_edge)/2-0.007,
                            height=0.31,
                            y=0.69,
                            vjust=1,
                            gp=gpar(lwd=1,fill=NA))
    dst_outline = rectGrob(width=1-text_edge-0.06,
                            x=text_edge + (1-text_edge)/2-0.007,
                            height=0.31,
                            y=0.18,
                            gp=gpar(lwd=1,fill=NA))

    western = gTree(children = gList(
        # alignment_lines,
        blot,
        antigen_labels,
        plus_labels,
        mutant_labels,
        temp_lines,
        temp_labels,
        spt6_outline,
        dst_outline))

    ggsave(pdf_out,
           plot=western,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     blot_path = snakemake@input[["blot_path"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

