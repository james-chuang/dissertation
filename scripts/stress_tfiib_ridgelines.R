
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

main = function(theme_spec,
                fonts_path,
                data_path,
                diamide_genic,
                diamide_intragenic,
                amino_genic,
                amino_intragenic,
                nitrogen_genic,
                nitrogen_intragenic,
                annotations,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)
    library(cowplot)
    library(ggridges)

    ttf_import(fonts_path)
    loadfonts()

    annotation_df = tibble()
    annotation_groups = c("DAN", "DA", "DN", "AN", "D", "A", "N")

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
        filter(group %in% c("YPD", "diamide", "SC", "SD", "nitrogen")) %>%
        mutate(annotation = fct_inorder(annotation, ordered=TRUE)) %>%
        group_by(group, annotation, index, position) %>%
        summarise(signal = mean(signal, na.rm=TRUE)) %>%
        # ungroup() %>%
        # complete(nesting(annotation, group, index),
        #          position,
        #          fill = list(signal = min(.[["signal"]])))
        group_by(annotation, index) %>%
        mutate(signal = scales::rescale(signal))

    peaks_lfc = build_peak_lfc_df(genic_results = diamide_genic,
                                  intragenic_results = diamide_intragenic,
                                  annotation_df = annotation_df,
                                  condition_id = "diamide") %>%
        bind_rows(build_peak_lfc_df(genic_results = amino_genic,
                                    intragenic_results = amino_intragenic,
                                    annotation_df = annotation_df,
                                    condition_id = "SD")) %>%
        bind_rows(build_peak_lfc_df(genic_results = nitrogen_genic,
                                    intragenic_results = nitrogen_intragenic,
                                    annotation_df = annotation_df,
                                    condition_id = "nitrogen"))

    df_lfc = df %>%
        left_join(peaks_lfc,
                  by = c("group", "annotation", "index")) %>%
        filter(position >= start_relative & position < end_relative) %>%
        select(-c(start_relative, end_relative))

    df %<>%
        left_join(df_lfc,
                  by = c("group", "annotation", "index", "position", "signal")) %>%
        mutate(group = ordered(group,
                               levels = c("YPD",
                                          "SC",
                                          "diamide",
                                          "SD",
                                          "nitrogen"),
                               labels = c("YPD",
                                          "SC",
                                          "oxidative stress",
                                          "amino acid stress",
                                          "nitrogen stress")))

    tfiib_heatmap = ggplot() +
        geom_ridgeline_gradient(data = df,
                                aes(x=position, y=index, height=5*signal,
                                    fill=log2_foldchange,
                                    group=index),
                       size=0.05) +
        facet_grid(annotation ~ group,
                   scales = "free_y",
                   space = "free_y") +
        scale_fill_distiller(palette = "PRGn",
                             limits = c(-2,2),
                             oob = scales::squish,
                             name = expression(log[2] ~ textstyle(frac(condition, control))),
                             labels = c(paste0("\U2264", "2"), "-1", "0", "1",
                                        paste0("\U2265", "2")),
                             guide=guide_colorbar(barwidth=10, barheight=0.8,
                                                  title.vjust=1),
                             na.value = "gray70") +
        scale_x_continuous(labels = function(x){case_when(x==0 ~ "ATG",
                                                          x==2 ~ "2 kb",
                                                          TRUE ~ as.character(x))},
                           expand = c(0, 0.05),
                           breaks = scales::pretty_breaks(n=3)) +
        scale_y_reverse(expand = c(0,0),
                        name = "TFIIB ChIP-nexus protection",
                        position = "right") +
        coord_cartesian(clip="off") +
        theme_heatmap +
        theme(axis.title.y = element_text(color="black", size=9, angle=90),
              panel.grid.major.x = element_line(color="gray80", size=0.1))
        # theme_heatmap +
        # theme(panel.border = element_rect(color="red", fill=NA))

    diagram_df = df %>%
        group_by(annotation) %>%
        summarise(index = median(index)-2.5) %>%
        mutate(D = if_else(str_detect(annotation, "D"), TRUE, FALSE),
               A = if_else(str_detect(annotation, "A"), TRUE, FALSE),
               N = if_else(str_detect(annotation, "N"), TRUE, FALSE)) %>%
        gather(key = "condition",
               value = "presence",
               -c(annotation, index)) %>%
        mutate(condition = ordered(condition,
                                   levels = c("D", "A", "N"),
                                   labels = c("oxidative stress",
                                              "amino acid stress",
                                              "nitrogen stress")))

    diagram = ggplot() +
        geom_hline(data = diagram_df,
                   aes(yintercept = index),
                   color="grey85",
                   size=0.5) +
        geom_point(data = diagram_df,
                   aes(x=condition,
                       y=index,
                       fill=presence),
                   shape=21,
                   size=3,
                   color="grey85") +
        geom_raster(data = df,
                   aes(x=1,
                       y=index),
                   fill="NA") +
        facet_grid(annotation~.,
                   scales = "free_y",
                   space = "free_y") +
        scale_fill_manual(values = c("grey95", "grey20"),
                          guide=FALSE) +
        scale_x_discrete(expand = c(0.2,0),
                         position = "top",
                         name = "induced in:") +
        scale_y_continuous(expand = c(0,1)) +
        coord_cartesian(clip="off") +
        theme_heatmap +
        theme(plot.margin = margin(r=-0.3, unit="cm"),
              axis.text.x.top = element_text(angle=15, vjust=0, hjust=0.1,
                                             margin = margin(b=5, unit="pt")),
              axis.title.x = element_text(size=9, color="black"))
              # axis.title.x = element_text(size=9, color="black"),
              # panel.border = element_rect(color="red", fill=NA))

    figure_1a = plot_grid(diagram, tfiib_heatmap, align = "h",
              axis = "tb", nrow=1,
              rel_widths = c(0.15, 1))

    ggplot2::ggsave(pdf_out, plot=figure_1a, width=fig_width, height=fig_height, units="in")
    embed_fonts(pdf_out)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path = snakemake@input[["fonts_path"]],
     data_path = snakemake@input[["data_path"]],
     annotations = snakemake@input[["annotations"]],
     diamide_genic = snakemake@input[["diamide_genic"]],
     diamide_intragenic = snakemake@input[["diamide_intragenic"]],
     amino_genic = snakemake@input[["amino_genic"]],
     amino_intragenic = snakemake@input[["amino_intragenic"]],
     nitrogen_genic = snakemake@input[["nitrogen_genic"]],
     nitrogen_intragenic = snakemake@input[["nitrogen_intragenic"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

