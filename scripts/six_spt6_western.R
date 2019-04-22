
main = function(theme_spec,
                fonts_path,
                data_path,
                blot_path,
                fig_width,
                fig_height,
                pdf_out){
    source(theme_spec)
    library(png)

    ttf_import(fonts_path)
    loadfonts()

    df = read_tsv(data_path) %>%
        select(1,2,3,4,5,18,19,20,21,22) %>%
        magrittr::set_colnames(c("strain", "temperature", "antigen", "background", "replicate",
                                 "gray_min", "grey_max", "gray_mean", "gray_median", "auc")) %>%
        select(strain, temperature, antigen, background, replicate, auc) %>%
        spread(background, auc) %>%
        mutate(signal = pmax(`FALSE`-`TRUE`, 0)) %>%
        select(-c(`FALSE`, `TRUE`)) %>%
        spread(antigen, signal) %>%
        mutate(spikenorm = Spt6/spikein) %>%
        group_by(strain, temperature) %>%
        mutate(group_mean = mean(spikenorm, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(strain = ordered(strain, levels = c("WT", "spt6-1004")))

    #rescale data so that mean of WT group is 1
    wt_og_mean = df %>%
        filter(strain=="WT", temperature==30) %>%
        distinct(group_mean) %>%
        pull(group_mean)

    df %<>%
        mutate(spikenorm_scaled = scales::rescale(spikenorm, from=c(0, wt_og_mean)))

    summary_df = df %>%
        group_by(strain, temperature) %>%
        summarise(group_mean_scaled = mean(spikenorm_scaled),
                  group_sd = sd(spikenorm_scaled))

    text_edge = 0.25
    right_edge=0.09
    increment = (1-text_edge-right_edge)/4
    start = increment/2+0.02

    alignment_lines = segmentsGrob(x0=seq(0,3)*increment+start+text_edge,
                                   x1=seq(0,3)*increment+start+text_edge,
                                   y0=0,
                                   y1=1,
                                   gp = gpar(lty="dashed"))
    temp_lines = segmentsGrob(x0=c(0,2)*increment-0.06+text_edge+start,
                                 x1=c(1,3)*increment+0.06+text_edge+start,
                                 y0=0.89, y1=0.89)
    temp_labels = textGrob(label = c(paste0("30", "\U00B0", "C"),
                                     paste0("37", "\U00B0", "C")),
                              gp = gpar(fontsize=9,
                                        fontfamily="FreeSans"),
                              x=text_edge+start+c(1,5)*increment/2,
                              y=1,
                              vjust=1.1)
    antigen_labels = textGrob(label=c("Spt6-FLAG",
                                      "Dst1-Myc",
                                      "Spt6:"),
                              x=text_edge,
                              y=c(0.61,
                                  0.36,
                                  0.82),
                              gp=gpar(fontsize=9,
                                      fontfamily="FreeSans"),
                              hjust=1)
    plus_labels = textGrob(label = "+",
                           x=text_edge+start+c(0,2)*increment,
                           y=0.82,
                           gp=gpar(fontsize=9,
                                   fontfamily="FreeSans"),
                           vjust=0.5)
    mutant_labels = textGrob(label = "1004",
                             x=text_edge+start+c(1,3)*increment,
                             y=0.82,
                             gp=gpar(fontsize=9,
                                     fontface="italic",
                                     fontfamily="FreeSans"),
                             vjust=0.5)
    blot = rasterGrob(readPNG(blot_path),
                      width=1-text_edge+0.04,
                      x=text_edge + (1-text_edge)/2,
                      height=0.65,
                      y=unit(0.5, "npc"))
    spt6_outline = rectGrob(width=1-text_edge-0.06,
                            x=text_edge + (1-text_edge)/2-0.007,
                            height=0.23,
                            y=unit(0.63, "npc"),
                            gp=gpar(lwd=1,fill=NA))
    dst_outline = rectGrob(width=1-text_edge-0.06,
                            x=text_edge + (1-text_edge)/2-0.007,
                            height=0.23,
                            y=unit(0.35, "npc"),
                            gp=gpar(lwd=1,fill=NA))
    quant_labels = textGrob(label = c(expression(textstyle(atop("1.00", phantom(.) %+-% 0.25 ~ phantom(.)))),
                                      expression(textstyle(atop(0.75, phantom(.) %+-% 0.12 ~ phantom(.)))),
                                      expression(textstyle(atop(0.84, phantom(.) %+-% 0.01 ~ phantom(.)))),
                                      expression(textstyle(atop(0.19, phantom(.) %+-% 0.05 ~ phantom(.))))),
                            gp=gpar(fontsize=9,
                                    fontfamily="FreeSans"),
                            x=text_edge+start+seq(0,3)*increment,
                            y=0,
                            vjust=0)

    western = gTree(children = gList(
        # alignment_lines,
        blot,
        antigen_labels,
        plus_labels,
        mutant_labels,
        temp_lines,
        temp_labels,
        spt6_outline, dst_outline,
        quant_labels))

    ggsave(pdf_out,
           plot=western,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path = snakemake@input[["fonts_path"]],
     data_path = snakemake@input[["data_path"]],
     blot_path = snakemake@input[["blot_path"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

