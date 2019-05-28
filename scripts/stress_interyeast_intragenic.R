
main = function(theme_spec = "thesis_theme.R",
                mik_input_path = "scer_intragenic_peaks-v-smik.tsv",
                uva_input_path = "scer_intragenic_peaks-v-sbay.tsv",
                pdf_out = "test.pdf",
                fig_height = 6,
                fig_width = 2) {

    source(theme_spec)

    df = read_tsv(mik_input_path) %>%
        transmute(scer_peak_id = tss_name,
                  scer_lfc = tss_lfc,
                  # scer_lfc_se = tss_lfc_SE,
                  scer_fdr = tss_FDR,
                  scer_feature_id = feature_name,
                  smik_peak_id = reference_name,
                  smik_lfc = as.numeric(reference_lfc),
                  smik_fdr = as.numeric(reference_qval)) %>%
        inner_join(read_tsv(uva_input_path)%>%
                       transmute(scer_peak_id = tss_name,
                                 scer_lfc = tss_lfc,
                                 # scer_lfc_se = tss_lfc_SE,
                                 scer_fdr = tss_FDR,
                                 sbay_peak_id = reference_name,
                                 sbay_lfc = as.numeric(reference_lfc),
                                 sbay_fdr = as.numeric(reference_qval)),
                   by = c("scer_peak_id",
                          "scer_lfc",
                          "scer_fdr")) %>%
        mutate(max_other_fdr = pmax(smik_fdr, sbay_fdr, na.rm=TRUE))

    heatmap_df = df %>%
        filter(scer_fdr > 1) %>%
        group_by(scer_peak_id, scer_lfc, scer_fdr, scer_feature_id) %>%
        arrange(desc(max_other_fdr)) %>%
        slice(1) %>%
        ungroup() %>%
        arrange(scer_lfc) %>%
        mutate(scer_peak_id = fct_inorder(scer_peak_id, ordered=TRUE)) %>%
        select(scer_peak_id,
               scer_feature_id,
               scer_lfc,
               smik_lfc,
               sbay_lfc) %>%
        gather(key = species,
               value = lfc,
               -c(scer_peak_id, scer_feature_id)) %>%
        distinct() %>%
        mutate(species = ordered(species,
                                 levels = c("sbay_lfc",
                                           "smik_lfc",
                                           "scer_lfc"),
                                 labels = c("S. bayanus",
                                           "S. mikatae",
                                           "S. cerevisiae")))

    plot = ggplot(data = heatmap_df,
                  aes(y = species,
                      x = scer_peak_id,
                      fill = lfc)) +
        geom_tile(interpolate=FALSE) +
        scale_fill_distiller(palette = "PRGn",
                             limits = c(-5, 5),
                             oob = scales::squish,
                             breaks = c(-5,-2.5,0,2.5,5),
                             labels = c("≤-5","-2.5","0","2.5","≥5"),
                             name = expression("log"[2] ~ textstyle(frac("oxidative stress", "YPD"))),
                             guide = guide_colorbar(barwidth=10,
                                                    title.vjust=1)) +
        scale_x_discrete(labels = heatmap_df[["scer_feature_id"]],
                         position="top",
                         expand=c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        ggtitle(expression(italic("S. cerevisiae") * " oxidative-stress-dependent intragenic TSSs")) +
        theme_default +
        theme(axis.text.x = element_text(angle=50,
                                         size=8,
                                         hjust=0,
                                         color="black"),
              axis.text.y = element_text(face="italic"),
              panel.grid.major.x = element_blank(),
              axis.ticks.x = element_line(color="black"),
              axis.title = element_blank(),
              legend.position="bottom",
              legend.title=element_text(size=10),
              legend.justification=c(0.5,0.5),
              plot.margin=margin(r=8, unit="pt"),
              panel.border=element_blank())

    ggsave(pdf_out,
           plot = plot,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     mik_input_path = snakemake@input[["mikatae"]],
     uva_input_path = snakemake@input[["uvarum"]],
     pdf_out = snakemake@output[["pdf"]],
     fig_height = snakemake@params[["height"]],
     fig_width = snakemake@params[["width"]])

