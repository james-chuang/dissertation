
import = function(path, group){
    read_tsv(path) %>%
        transmute(occupancy = smt_value,
                  fuzziness = fuzziness_score,
                  group = group) %>%
        return()
}


main = function(theme_spec = "thesis_theme.R",
                fonts_path,
                wt_mnase_quant = "nucleosome_quantification_data_non-depleted.Fnor.smooth.positions.xls",
                mut_mnase_quant = "nucleosome_quantification_data_depleted.Fnor.smooth.positions.xls",
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)

    ttf_import(fonts_path)
    loadfonts()

    df = import(wt_mnase_quant, group="non-depleted") %>%
        bind_rows(import(mut_mnase_quant, group="depleted")) %>%
        mutate(group = fct_inorder(group, ordered=TRUE))

    summary_df = df %>%
        group_by(group) %>%
        summarise(mean_occ = mean(occupancy),
                  sd_occ = sd(occupancy),
                  mean_fuzz = mean(fuzziness),
                  sd_fuzz = sd(fuzziness),
                  median_occ = median(occupancy),
                  median_fuzz = median(fuzziness))

    figure = ggplot() +
        geom_density(data = df,
                     aes(x=fuzziness,
                         fill=group,
                         color=group),
                     alpha=0.6,
                     bw = 0.7) +
        scale_x_continuous(breaks = scales::pretty_breaks(n=3),
                           expand = c(0,0),
                           name = expression("fuzziness" %==% "std. dev of dyad positions (bp)")) +
        scale_y_continuous(breaks = scales::pretty_breaks(n=2),
                           expand = c(0,0)) +
        scale_fill_few() +
        scale_color_few() +
        theme_default +
        theme(panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank()) +
        theme(legend.position = c(0.33, 0.5))

    ggsave(pdf_out,
           plot=figure,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path = snakemake@input[["fonts_path"]],
     wt_mnase_quant = snakemake@input[["wt_mnase_quant"]],
     mut_mnase_quant = snakemake@input[["mut_mnase_quant"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

