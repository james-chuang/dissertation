
main = function(theme_spec = "thesis_theme.R",
                dsk_scer = "scer_dsk2.png",
                snr_scer = "scer_u190.png",
                hsp_scer = "scer_hsp12.png",
                dsk_sbay = "sbay_dsk2.png",
                snr_sbay = "sbay_u190.png",
                hsp_sbay = "sbay_hsp12.png",
                dsk_smik = "smik_dsk2.png",
                snr_smik = "smik_u190.png",
                hsp_smik = "smik_hsp12.png",
                fig_width = 5,
                fig_height = 1.5,
                pdf_out){
    source(theme_spec)
    library(png)

    left_margin = 0.11
    right_margin= 0.23
    blot_buffer = 0.01

    blot_width = (1-right_margin-left_margin-2*blot_buffer)/3

    species_labels = textGrob(label=c("S. cerevisiae",
                                      "S. mikatae",
                                      "S. bayanus"),
                               x=c(left_margin+blot_width/2,
                                   left_margin+3*blot_width/2+blot_buffer,
                                   left_margin+5*blot_width/2+2*blot_buffer),
                               y=1,
                               gp=gpar(fontsize=10,
                                       fontface="italic",
                                       fontfamily="FreeSans"),
                               hjust=0.5,
                               vjust=1)

    species_lines = segmentsGrob(x0=c(left_margin + 0.10*blot_width,
                                      left_margin+blot_width+blot_buffer + 0.10*blot_width,
                                      left_margin+2*(blot_width+blot_buffer) + 0.10*blot_width),
                                  x1=c(left_margin+blot_width - 0.10*blot_width,
                                       left_margin+blot_buffer+2*blot_width - 0.10*blot_width,
                                       left_margin+2*blot_buffer+3*blot_width - 0.10*blot_width),
                                  y0=unit(1, "npc") - unit(1.3, "strheight", data="S"),
                                  y1=unit(1, "npc") - unit(1.3, "strheight", data="S"))

    diamide_labels = textGrob(label=c("diamide:", "-", "+", "-", "+", "-", "+"),
                              x=c(left_margin/2,
                                  left_margin+blot_width/2 - 0.23*blot_width,
                                  left_margin+blot_width/2 + 0.23*blot_width,
                                  left_margin+3*blot_width/2+blot_buffer - 0.23*blot_width,
                                  left_margin+3*blot_width/2+blot_buffer + 0.23*blot_width,
                                  left_margin+5*blot_width/2+2*blot_buffer - 0.23*blot_width,
                                  left_margin+5*blot_width/2+2*blot_buffer + 0.23*blot_width),
                              y=unit(1, "npc") - unit(1.8, "strheight", data="S"),
                              gp=gpar(fontsize=10,
                                      fontfamily="FreeSans"),
                              hjust=0.5,
                              vjust=0.5)

    transcript_labels = textGrob(label=c(expression(italic("DSK2") * " full-length"),
                                         expression(italic("DSK2") * " intragenic"),
                                         expression(italic("HSP12")),
                                         expression(italic("SNR190"))),
                                 x=1-right_margin*0.95,
                                 y=c(0.70,
                                     0.54,
                                     0.32,
                                     0.10),
                                 gp=gpar(fontsize=10,
                                         fontfamily="FreeSans"),
                                 hjust=0)

    dsk_scer_blot = rasterGrob(readPNG(dsk_scer)[20:265,
                                                 40:270],
                               width=blot_width,
                               height = (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.5,
                               x=left_margin+blot_width/2,
                               y=unit(1, "npc") - unit(2.3, "strheight", data="S"),
                               vjust=1)
    dsk_smik_blot = rasterGrob(readPNG(dsk_smik)[70:315,
                                                 65:295],
                             width=blot_width,
                             height = (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.5,
                             x=left_margin+blot_width*3/2+blot_buffer,
                             y=unit(1, "npc") - unit(2.3, "strheight", data="S"),
                             vjust=1)
    dsk_sbay_blot = rasterGrob(readPNG(dsk_sbay)[50:295,
                                                 30:260],
                             width=blot_width,
                             height = (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.5,
                             x=left_margin+blot_width*5/2+2*blot_buffer,
                             y=unit(1, "npc") - unit(2.3, "strheight", data="S"),
                             vjust=1)
    hsp_scer_blot = rasterGrob(readPNG(hsp_scer)[10:100,
                                                 20:260],
                               width=blot_width,
                               height = (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.25,
                               x=left_margin+blot_width/2,
                               y=(unit(1, "npc") - unit(2.3, "strheight", data="S")) -
                                      (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.5 -
                                   unit(blot_buffer, "npc")*3,
                               vjust=1)
    hsp_smik_blot = rasterGrob(readPNG(hsp_smik)[30:120,
                                                 30:270],
                             width=blot_width,
                             height = (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.25,
                             x=left_margin+blot_width*3/2+blot_buffer,
                             y=(unit(1, "npc") - unit(2.3, "strheight", data="S")) -
                                 (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.5 -
                                 unit(blot_buffer, "npc")*3,
                             vjust=1)
    hsp_sbay_blot = rasterGrob(readPNG(hsp_sbay)[10:100,
                                                 40:280],
                             width=blot_width,
                             height = (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.25,
                             x=left_margin+blot_width*5/2+2*blot_buffer,
                             y=(unit(1, "npc") - unit(2.3, "strheight", data="S")) -
                                 (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.5 -
                                 unit(blot_buffer, "npc")*3,
                             vjust=1)
    snr_scer_blot = rasterGrob(readPNG(snr_scer)[10:140,
                                                 25:265],
                               width=blot_width,
                               height = (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.25,
                               x=left_margin+blot_width/2,
                               y=(unit(1, "npc") - unit(2.3, "strheight", data="S")) -
                                      (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.75 -
                                   unit(blot_buffer, "npc")*6,
                               vjust=1)
    snr_smik_blot = rasterGrob(readPNG(snr_smik)[35:165,
                                                 45:285],
                             width=blot_width,
                             height = (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.25,
                             x=left_margin+blot_width*3/2+blot_buffer,
                             y=(unit(1, "npc") - unit(2.3, "strheight", data="S")) -
                                 (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.75 -
                                 unit(blot_buffer, "npc")*6,
                             vjust=1)
    snr_sbay_blot = rasterGrob(readPNG(snr_sbay)[5:135,
                                                 25:265],
                             width=blot_width,
                             height = (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.25,
                             x=left_margin+blot_width*5/2+2*blot_buffer,
                             y=(unit(1, "npc") - unit(2.3, "strheight", data="S")) -
                                 (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc"))*0.75 -
                                 unit(blot_buffer, "npc")*6,
                             vjust=1)

    blot_outlines = rectGrob(width=3*blot_width+2*blot_buffer,
                             x=left_margin+3/2*blot_width+blot_buffer,
                             height = (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc")) * c(0.5,0.25,0.25),
                             y = (unit(1, "npc") - unit(2.3, "strheight", data="S")) -
                                 (unit(1, "npc") - unit(2.3, "strheight", data="S") - unit(7*blot_buffer, "npc")) *
                                      c(0, 0.5, 0.75) -
                                 unit(blot_buffer, "npc") *
                                 c(0, 3, 6),
                             vjust=1,
                             gp=gpar(lwd=1,fill=NA))


    northern = gTree(children = gList(species_labels,
                                      species_lines,
                                      diamide_labels,
                                      transcript_labels,
                                      dsk_scer_blot,
                                      dsk_smik_blot,
                                      dsk_sbay_blot,
                                      hsp_scer_blot,
                                      hsp_smik_blot,
                                      hsp_sbay_blot,
                                      snr_scer_blot,
                                      snr_smik_blot,
                                      snr_sbay_blot,
                                      blot_outlines
                                      ))
    ggsave(pdf_out,
           plot=northern,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)

}

main(theme_spec = snakemake@input[["theme"]],
     dsk_scer = snakemake@input[["dsk_scer"]],
     snr_scer = snakemake@input[["snr_scer"]],
     hsp_scer = snakemake@input[["hsp_scer"]],
     dsk_smik = snakemake@input[["dsk_smik"]],
     snr_smik = snakemake@input[["snr_smik"]],
     hsp_smik = snakemake@input[["hsp_smik"]],
     dsk_sbay = snakemake@input[["dsk_sbay"]],
     snr_sbay = snakemake@input[["snr_sbay"]],
     hsp_sbay = snakemake@input[["hsp_sbay"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

