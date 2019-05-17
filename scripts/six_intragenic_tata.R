
import = function(path, annotation_id){
    read_tsv(path, col_types = 'ciicdcciicdccc') %>%
        mutate(start = case_when(region_strand=="+" & motif_chrom != "." ~ motif_start-region_end,
                                 region_strand=="-" & motif_chrom != "." ~ region_start-motif_end,
                                 motif_chrom=="." ~ as.integer(NaN)),
               end = case_when(region_strand=="+" & motif_chrom != "." ~ motif_end-region_end,
                               region_strand=="-" & motif_chrom != "." ~ region_start-motif_start,
                               motif_chrom=="." ~ as.integer(NaN)),
               annotation = annotation_id) %>%
        arrange(region_score) %>%
        return()
}

main = function(theme_spec,
                tata_genic_path,
                tata_intra_path,
                tata_anti_path,
                # tata_random_path,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)
    library(broom)

    df = tata_anti_path %>%
        import(annotation_id = "antisense") %>%
        bind_rows(tata_genic_path %>% import(annotation_id = "genic")) %>%
        bind_rows(tata_intra_path %>% import(annotation_id = "intragenic")) %>%
        group_by(annotation) %>%
        mutate(n_regions = n_distinct(chrom, region_start, region_end, region_strand)) %>%
        filter(! is.na(start)) %>%
        mutate(n_motifs = n()) %>%
        ungroup() %>%
        transmute(annotation = fct_inorder(annotation, ordered=TRUE),
                  n_regions = n_regions,
                  n_motifs = n_motifs,
                  position = (end+start)/2) %>%
        group_by(annotation, n_regions, n_motifs) %>%
        do(density(x=.$position,
                   kernel="gaussian",
                   bw=3.333) %>%
               tidy())

    tata = ggplot(data = df,
           aes(x=x,
               y=y*n_motifs/n_regions,
               group=annotation,
               fill=annotation)) +
        geom_area(na.rm=TRUE,
                  # fill=NA
                  alpha=0.8,
                  linetype="blank",
                  size=0) +
        geom_vline(xintercept = 0, color="grey65", size=1.5) +
        scale_x_continuous(limits = c(-200, 0),
                           expand = c(0,0),
                           breaks = scales::pretty_breaks(n=3),
                           labels = function(x) case_when(x==-200 ~ "-200 nt",
                                                          x==0 ~ "TSS",
                                                          TRUE ~ as.character(x))) +
        scale_y_continuous(expand = c(0,0),
                           breaks = c(0,0.008),
                           name = "scaled density") +
        scale_fill_manual(values = c("#E15759",
                                     "#4E79A7",
                                     "#F28E2B")) +
        ggtitle("TATA consensus probability") +
        theme_default +
        theme(axis.title.x = element_blank(),
              axis.line.y = element_line(size=0.25, color="grey65"),
              axis.text.x = element_text(size=10),
              axis.text.y = element_text(size=8),
              axis.title.y = element_text(margin = margin(r=-5, unit="pt")),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              legend.justification = c(0.5, 0.5),
              legend.position = c(0.30, 0.55),
              legend.key.width= unit(12, "pt"),
              legend.spacing.x = unit(1, "pt"),
              plot.margin = margin(0,10,0,0,"pt"))

    ggsave(pdf_out,
           plot=tata,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     tata_genic_path = snakemake@input[["tata_genic_path"]],
     tata_intra_path = snakemake@input[["tata_intra_path"]],
     tata_anti_path = snakemake@input[["tata_anti_path"]],
     # tata_random_path = snakemake@input[["tata_random_path"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])
