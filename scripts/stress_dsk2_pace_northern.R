
main = function(theme_spec = "thesis_theme.R",
                dsk_blot = "JLW_20190509 DSK2-[Phosphor]_rotated.jpg",
                fig_width,
                fig_height,
                pdf_out){
    source(theme_spec)
    library(jpeg)

    left_margin = 0.17
    right_margin = 0.355
    blot_buffer = 0.01

    blot_width = (1-right_margin-left_margin-blot_buffer)/2

    dsk_image = readJPEG(dsk_blot)

    wildtype_blot = rasterGrob(dsk_image[1450:2000,
                                         750:1230,],
                               width=blot_width,
                               height=0.56,
                               x=left_margin+blot_width/2,
                               y=0.25,
                               vjust=0)
    mutant_blot = rasterGrob(dsk_image[1450:2000,
                                       1620:2100,],
                             width=blot_width,
                             height=0.56,
                             x=left_margin+blot_width*3/2+blot_buffer,
                             y=0.25,
                             vjust=0)
    blot_outlines = rectGrob(width=2*blot_width+blot_buffer,
                             x=left_margin+blot_width+blot_buffer/2,
                             height=c(0.56, 0.21),
                             y=c(0.25,0.01),
                             vjust=0,
                             gp=gpar(lwd=1,fill=NA))

    genotype_labels = textGrob(label=c("DSK2",
                                       "dsk2-pace"),
                               x=c(left_margin+blot_width/2,
                                   left_margin+3*blot_width/2+blot_buffer),
                               y=1,
                               gp=gpar(fontsize=10,
                                       fontface="italic",
                                       fontfamily="FreeSans"),
                               hjust=0.5,
                               vjust=1)

    genotype_lines = segmentsGrob(x0=c(left_margin + 0.10*blot_width,
                                       left_margin+blot_width+blot_buffer + 0.10*blot_width),
                                  x1=c(left_margin+blot_width - 0.10*blot_width,
                                       left_margin+blot_buffer+2*blot_width - 0.10*blot_width),
                                  y0=unit(1, "npc") - unit(1.3, "strheight", data="DSK2"),
                                  y1=unit(1, "npc") - unit(1.3, "strheight", data="DSK2"))

    diamide_labels = textGrob(label=c("diamide:", "-", "+", "-", "+"),
                              x=c(left_margin,
                                  left_margin+blot_width/2 - 0.23*blot_width,
                                  left_margin+blot_width/2 + 0.23*blot_width,
                                  left_margin+3*blot_width/2+blot_buffer - 0.23*blot_width,
                                  left_margin+3*blot_width/2+blot_buffer + 0.23*blot_width),
                              y=unit(1, "npc") - unit(1.8, "strheight", data="DSK2"),
                              gp=gpar(fontsize=10,
                                      fontfamily="FreeSans"),
                              hjust=1,
                              vjust=0.5)
    transcript_labels = textGrob(label=c(expression(italic("DSK2") * " full-length"),
                                         expression(italic("DSK2") * " intragenic"),
                                         "control"),
                                 x=1-right_margin*0.95,
                                 y=c(0.66,
                                     0.38,
                                     0.115),
                                 gp=gpar(fontsize=10,
                                         fontfamily="FreeSans"),
                                 hjust=0)

    temp_label = textGrob(label="loading control coming",
                          x=left_margin+blot_width+blot_buffer/2,
                          y=0.115,
                          gp=gpar(fontsize=10,
                                  fontfamily="FreeSans")
                          )

    # alignment_line = segmentsGrob(x0=0,
    #                               x1=1,
    #                               y0=0.115,
    #                               y1=0.115)

    northern = gTree(children = gList(wildtype_blot,
                                      mutant_blot,
                                      genotype_lines,
                                      genotype_labels,
                                      diamide_labels,
                                      transcript_labels,
                                      # alignment_line,
                                      temp_label,
                                      blot_outlines
                                      ))
    # grid.newpage()
    # grid.draw(northern)

    ggsave(pdf_out,
           plot=northern,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     dsk_blot = snakemake@input[["dsk_blot"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

