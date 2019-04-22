
import = function(df,
                  path,
                  label_col_id,
                  category){
    df = read_tsv(path) %>%
        distinct(chrom,
                 start,
                 end,
                 name,
                 score,
                 strand,
                 .keep_all = TRUE) %>%
        select(chrom,
               start,
               end,
               label=label_col_id,
               score,
               strand,
               log2_foldchange,
               lfc_SE,
               stat,
               log10_pval,
               log10_padj,
               mean_expr,
               condition_expr,
               control_expr) %>%
        mutate(category = category) %>%
        bind_rows(df, .) %>%
        return()
}

main = function(theme_spec,
                fonts_path,
                in_genic,
                in_intra,
                in_anti,
                in_inter,
                alpha,
                fig_width,
                fig_height,
                pdf_out){
    source(theme_spec)

    ttf_import(fonts_path)
    loadfonts()

    df = tibble() %>%
        import(in_genic, label_col_id="genic_name", category="genic") %>%
        import(in_intra, label_col_id="orf_name", category="intragenic") %>%
        import(in_anti, label_col_id="transcript_name", category="antisense") %>%
        import(in_inter, label_col_id="name", category="intergenic") %>%
        mutate(category = fct_inorder(category, ordered=TRUE),
               significance = if_else(log10_padj > -log10(alpha), TRUE, FALSE)) %>%
        mutate(direction = if_else(log2_foldchange >= 0, "upregulated", "downregulated"))

    summary_df = df %>% count(category, significance, direction) %>%
        group_by(category) %>%
        mutate(group_total = sum(n)) %>%
        ungroup() %>%
        spread(direction, n) %>%
        filter(significance) %>%
        mutate(ymax=cumsum(group_total),
               ymin=ymax-group_total,
               y=(ymax+ymin)/2) %>%
        rowid_to_column(var = "ymax_unscaled") %>%
        mutate(ymin_unscaled = lag(ymax_unscaled, default=as.integer(0)),
               y_unscaled = (ymin_unscaled+ymax_unscaled)/2)

    diffexp_summary = ggplot(data = summary_df) +
        geom_segment(aes(x=if_else(category=="genic", upregulated, -downregulated),
                         # xend=if_else(category=="genic", 1700, -3000),
                         xend=-300,
                         y=y+500, yend=y+500),
                         # y=y_unscaled+0.17, yend=y_unscaled+0.17),
                         alpha=0.2,
                         size=0.2) +
        geom_rect(aes(xmin=0, xmax=upregulated,
                      ymax=ymax, ymin=ymin),
                      # ymax=ymax_unscaled, ymin=ymin_unscaled),
                  # size=0.5, color="white", fill=ptol_pal()(2)[2], alpha=0.7) +
                  size=0.5, color="white", fill="#af8dc3", alpha=0.9) +
        geom_rect(aes(xmin=-downregulated, xmax=0,
                      ymax=ymax, ymin=ymin),
                      # ymax=ymax_unscaled, ymin=ymin_unscaled),
                  size=0.5, color="white", fill="#7fbf7b", alpha=0.9) +
        geom_text(aes(x=upregulated/2,
                      y=y,
                      # y=y_unscaled,
                      label=upregulated),
                  hjust=0.5,
                  size=9/72*25.4,
                  family="FreeSans") +
        geom_text(aes(x=if_else(downregulated>20,
                                -downregulated/2,
                                -(downregulated+10)),
                      hjust=if_else(downregulated>20,
                                    0.5,
                                    1),
                      y=y,
                      # y=y_unscaled,
                      label=downregulated),
                  size=9/72*25.4,
                  family="FreeSans") +
        geom_text(aes(y=y,
                      # y=y_unscaled,
                      label=category),
                  x=-300,
                  hjust=0,
                  size=7/72*25.4,
                  family="FreeSans") +
        annotate(geom="label",
                 x=max(summary_df[["upregulated"]])/2,
                 fill= "#af8dc3",
                 label.r = unit(0, "pt"),
                 label.size = NA,
                 y=-1000,
                 # y=-0.3,
                 hjust=0.5,
                 label="upregulated",
                 size=7/72*25.4,
                 alpha=0.7,
                 family="FreeSans") +
        annotate(geom="label",
                 x=max(summary_df[["downregulated"]])/-2,
                 y=-1000,
                 # y=-0.3,
                 fill= "#7fbf7b",
                 label.r = unit(0, "pt"),
                 label.size = NA,
                 hjust=0.5,
                 label="downregulated",
                 size=7/72*25.4,
                 alpha=0.7,
                 family="FreeSans") +
        scale_y_reverse(limits = c(max(summary_df[["ymax"]]), -2000),
                        expand = c(0,0)) +
        # scale_y_reverse(limits = c(max(summary_df[["ymax_unscaled"]]), -.5),
        #                 expand = c(0,0))+
        theme_void() +
        theme(plot.margin = margin(0,0,0,0,"pt"))

    ggsave(pdf_out, plot=diffexp_summary, width=fig_width, height=fig_height, units="in")
    embed_fonts(pdf_out)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path = snakemake@input[["fonts_path"]],
     in_genic = snakemake@input[["in_genic"]],
     in_intra = snakemake@input[["in_intra"]],
     in_anti = snakemake@input[["in_anti"]],
     in_inter = snakemake@input[["in_inter"]],
     alpha = snakemake@params[["fdr_cutoff_tss"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

