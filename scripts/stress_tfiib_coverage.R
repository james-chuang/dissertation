import = function(df,
                  path,
                  gene_id){
    read_tsv(path,
             col_names = c("group", "sample",
                           "annotation", "assay",
                           "index", "position",
                           "signal")) %>%
        select(group, sample, position, signal) %>%
        mutate(gene=gene_id) %>%
        group_by(group, position, gene) %>%
        summarise(signal=mean(signal)) %>%
        bind_rows(df, .) %>%
        return()
}

import_bed = function(path){
    read_tsv(path,
             col_names = c("chrom", "start", "end",
                           "name", "score", "strand")) %>%
        select(-c(chrom, score)) %>%
        return()
}

main = function(theme_spec,
                fonts_path,
                data_paths,
                transcript_annotation, orf_annotation,
                gene_ids,
                up_distances, down_distances,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)
    library(cowplot)

    ttf_import(fonts_path)
    loadfonts()

    annotation = import_bed(transcript_annotation) %>%
        left_join(import_bed(orf_annotation),
                  by="name",
                  suffix=c("_transcript", "_orf")) %>%
        filter(name %in% gene_ids) %>%
        transmute(gene_id = ordered(name, levels=gene_ids),
                  orf_start = if_else(strand_transcript=="+",
                                      start_orf - start_transcript,
                                      end_transcript - end_orf),
                  orf_end = if_else(strand_transcript=="+",
                                    end_orf - start_transcript,
                                    end_transcript - start_orf),
                  transcript_end = abs(end_transcript - start_transcript)) %>%
        arrange(gene_id) %>%
        mutate(data_start = -up_distances,
               data_end = transcript_end + down_distances)

    gene_diagrams = ggplot(data = annotation) +
        geom_segment(aes(x=data_start, xend=data_end),
                     y=0, yend=0,
                     size=0, color=NA, fill=NA) +
        geom_segment(aes(xend=transcript_end),
                     x=0, y=0, yend=0, color="grey20") +
        geom_polygon(data = annotation %>%
                         select(-c(data_start, data_end, transcript_end)) %>%
                         mutate(notch = orf_end*0.93) %>%
                         gather(key, x, -gene_id) %>%
                         mutate(y = case_when(key=="orf_start" ~ list(c(-1,1)),
                                              key=="orf_end" ~ list(0),
                                              key=="notch" ~ list(c(1,-1)))) %>%
                         unnest(y) %>%
                         arrange(gene_id, key, x, y) %>%
                         group_by(gene_id) %>%
                         mutate(point_order = c(4,2,3,5,1)) %>%
                         arrange(point_order),
                         aes(x=x, y=y),
                     fill="grey80") +
        geom_text(aes(x=(orf_start+orf_end)/2,
                      label=paste0("italic(\"", gene_id, "\")")),
                  y=0, size=10/72*25.4, parse=TRUE,
                  family="FreeSans") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(limits = c(-1.5, 1.5), expand=c(0,0)) +
        facet_grid(.~gene_id,
                   scales="free_x",
                   space="free_x") +
        theme_void() +
        theme(strip.text.x = element_blank(),
              plot.margin = margin(b=-2, unit="pt"))

    df = tibble()
    for (i in 1:length(data_paths)){
        df %<>% import(data_paths[i], gene_ids[ceiling(i/3)])
    }
    df %<>%
        distinct() %>%
        mutate(group = ordered(group,
                               levels = c("untreated",
                                          "SC",
                                          "diamide",
                                          "SD",
                                          "nitrogen"),
                               labels = c("YPD",
                                          "SC",
                                          "oxidative\nstress",
                                          "amino acid\nstress",
                                          "nitrogen\nstress")),
               gene = fct_inorder(gene)) %>%
        group_by(group, gene, position) %>%
        summarise(signal = mean(signal)) %>%
        group_by(gene) %>%
        mutate(signal = scales::rescale(signal))

    coverage_plot = ggplot(data = df,
           aes(x=position, y=signal)) +
        geom_area(size=0.3, fill="#114477", color="#114477", alpha=0.9) +
        facet_grid(group ~ gene,
                   scales="free_x",
                   space="free_x",
                   switch="y") +
        scale_x_continuous(name = NULL,
                           labels = function(x) case_when(x==0 ~ "TSS",
                                                          x %in% c(1,2) ~ paste(x, "kb"),
                                                          TRUE ~ as.character(x)),
                           breaks = scales::pretty_breaks(n=2),
                           expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0),
                           breaks = scales::pretty_breaks(n=1),
                           sec.axis = dup_axis(name = "relative signal")) +
        theme_default +
        theme(strip.placement = "outside",
              strip.text = element_text(),
              strip.text.x = element_blank(),
              strip.text.y = element_text(angle=-180, hjust=1, vjust=0.5),
              axis.title.x = element_blank(),
              axis.title.y.left = element_blank(),
              axis.line.y.right = element_blank(),
              axis.ticks.y.right = element_blank(),
              axis.text.y = element_text(size=8),
              axis.text.y.right = element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_line(color="grey93"),
              panel.grid.major.y = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line.y = element_line(size=0.2, color="grey70"),
              plot.margin = margin(t=-2, unit="pt"))

    figure_1b = plot_grid(gene_diagrams, coverage_plot,
              align = "v", axis="lr",
              ncol=1,
              rel_heights = c(0.1, 1))

    ggplot2::ggsave(pdf_out,
                    plot=figure_1b,
                    width=fig_width,
                    height=fig_height,
                    units="in",
                    device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path = snakemake@input[["fonts_path"]],
     data_paths = snakemake@input[["data_paths"]],
     transcript_annotation = snakemake@input[["transcript_annotation"]],
     orf_annotation = snakemake@input[["orf_annotation"]],
     gene_ids = snakemake@params[["gene_ids"]],
     up_distances = snakemake@params[["upstream"]] ,
     down_distances = snakemake@params[["downstream"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

