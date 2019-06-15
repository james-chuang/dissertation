
import_annotation = function(df, input_path, annotation_id){
    df %>%
        bind_rows(read_tsv(input_path,
                           col_names = c("chrom",
                                         "start",
                                         "end",
                                         "name",
                                         "score",
                                         "strand")) %>%
                      rowid_to_column(var = "index") %>%
                      mutate(stop_position = (end-start)/1000,
                             annotation = annotation_id)) %>%
        return()
}

build_peak_lfc_df = function(genic_results,
                             intragenic_results,
                             annotation_df,
                             condition_id){
    read_tsv(genic_results) %>%
        select(start, end, log2_foldchange, gene_name = genic_name) %>%
        filter(gene_name %in% annotation_df[["gene_name"]]) %>%
        bind_rows(read_tsv(intragenic_results) %>%
                      select(start, end, log2_foldchange, gene_name = orf_name) %>%
                      filter(gene_name %in% annotation_df[["gene_name"]])) %>%
        left_join(annotation_df %>%
                      select(index, start, end, gene_name, strand, annotation),
                  by="gene_name",
                  suffix = c("", "_orf")) %>%
        mutate(start_relative = if_else(strand=="+", start-start_orf, end_orf-end),
               end_relative = if_else(strand=="+", end-start_orf, end_orf-start),
               group = condition_id) %>%
        mutate_at(vars(start_relative, end_relative),
                  funs(./1000)) %>%
        select(log2_foldchange, index, annotation, group, start_relative, end_relative) %>%
        return()
}

main = function(theme_spec = "thesis_theme.R",
                data_path = "genes-with-diamide-induced-intra-peaks-allsamples-allannotations-tfiib-chipnexus-libsizenorm-protection.tsv.gz",
                diamide_genic = "diamide-v-YPD_tfiib-chipnexus-libsizenorm-peaks-diffbind-results-genic-all.tsv",
                diamide_intragenic = "diamide-v-YPD_tfiib-chipnexus-libsizenorm-peaks-diffbind-results-intragenic-all.tsv",
                annotations = "orfs_w_diamide_induced_intragenic_TFIIB_peaks_ATGdistsort.bed",
                fig_width = 14,
                fig_height = 7.65,
                pdf_out = "test.pdf"){
    source(theme_spec)
    library(cowplot)
    library(ggridges)

    annotation_df = tibble()
    annotation_groups = c("D")

    for (annotation_index in 1:length(annotations)){
        annotation_df %<>%
            import_annotation(annotations[annotation_index],
                              annotation_groups[annotation_index])
    }
    annotation_df %<>%
        mutate(annotation = fct_inorder(annotation, ordered=TRUE)) %>%
        separate(name, into=c("gene_name", "peak_name"), extra="merge")

    df = read_tsv(data_path,
                  col_names = c("group",
                                "sample",
                                "annotation",
                                "index",
                                "position",
                                "signal")) %>%
        filter(group %in% c("YPD", "diamide")) %>%
        mutate(annotation = fct_inorder(annotation, ordered=TRUE)) %>%
        group_by(group, annotation, index, position) %>%
        summarise(signal = mean(signal, na.rm=TRUE)) %>%
        group_by(annotation, index) %>%
        mutate(signal = scales::rescale(signal))

    peaks_lfc = build_peak_lfc_df(genic_results = diamide_genic,
                                  intragenic_results = diamide_intragenic,
                                  annotation_df = annotation_df,
                                  condition_id = "diamide")

    df_lfc = df %>%
        left_join(peaks_lfc,
                  by = c("group", "index")) %>%
        filter(position >= start_relative & position < end_relative) %>%
        ungroup() %>%
        select(-c(start_relative, end_relative, annotation.x, annotation.y))

    df %<>%
        ungroup() %>%
        select(-annotation) %>%
        left_join(df_lfc,
                  by = c("group", "index", "position", "signal")) %>%
        mutate(group = ordered(group,
                               levels = c("YPD",
                                          "diamide"),
                               labels = c("unstressed",
                                          "oxidative stress")))

    tfiib_heatmap = ggplot() +
        geom_ridgeline_gradient(data = df,
                                aes(x=position, y=index, height=5*signal,
                                    fill=log2_foldchange,
                                    group=index),
                       size=0.05) +
        facet_grid(. ~ group,
                   scales = "free_y",
                   space = "free_y") +
        scale_fill_distiller(palette = "PRGn",
                             limits = c(-2,2),
                             oob = scales::squish,
                             name = expression("log"[2] * textstyle(frac("oxidative stress", "unstressed"))),
                             labels = c("≤-2", "-1", "0", "1","≥2"),
                             guide=guide_colorbar(barwidth=10, barheight=0.6,
                                                  title.vjust=1),
                             na.value = "gray70") +
        scale_x_continuous(labels = function(x){case_when(x==0 ~ "ATG",
                                                          x==2 ~ "2 kb",
                                                          TRUE ~ as.character(x))},
                           expand = c(0, 0.05),
                           breaks = scales::pretty_breaks(n=3)) +
        scale_y_reverse(expand = c(0,0),
                        name = paste(max(df[["index"]]), "induced intragenic peaks")) +
        ggtitle("relative TFIIB protection") +
        theme_heatmap_presentation +
        theme(panel.grid.major.x = element_line(color="gray80", size=0.1),
              legend.position="bottom",
              legend.text = element_text(size=10),
              legend.margin=margin(t=4, unit="pt"),
              strip.text = element_text(size=12,
                                        color="black",
                                        hjust=0.5,
                                        vjust=0.5,
                                        family="FreeSans"),
              strip.text.y = element_blank(),
              plot.title = element_text(size=12))

    ggplot2::ggsave(pdf_out,
                    plot=tfiib_heatmap,
                    width=fig_width,
                    height=fig_height,
                    units="cm",
                    device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     data_path = snakemake@input[["data_path"]],
     annotations = snakemake@input[["annotations"]],
     diamide_genic = snakemake@input[["diamide_genic"]],
     diamide_intragenic = snakemake@input[["diamide_intragenic"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])
