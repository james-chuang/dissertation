import_bed = function(path, species_id, annotation_type){
    read_tsv(path,
             col_names = c("chrom", "start", "end",
                           "name", "score", "strand")) %>%
        select(-c(chrom, score)) %>%
        mutate(species=species_id,
               annotation=annotation_type) %>%
        return()
}

main = function(theme_spec,
                transcripts_scer = "Scer_transcripts_w_verifiedORFs.bed",
                transcripts_smik = "Smik_codingtranscripts-adjustedTSS.bed",
                transcripts_sbay = "Sbay_codingtranscripts-adjustedTSS.bed",
                orfs_scer = "Scer_nondubious_ORFs_and_blocked_reading_frames-adjustedATG.bed",
                orfs_smik = "Smik_ORFs_3-adjustedATG.bed",
                orfs_sbay = "Sbay_ORFs_3-adjustedATG.bed",
                pace_scer = "DSK2_PACE_core_sites.bed",
                pace_smik = "Smik_DSK2_PACE_core_sites.bed",
                pace_sbay = "Sbay_DSK2_PACE_core_sites.bed",
                peak_scer = "diamide-v-YPD_tss-seq-spikenorm-peaks-diffexp-results-intragenic-all.tsv",
                peak_smik = "Scer-diamide-induced-TSS-Smik-coords.tsv",
                peak_sbay = "Scer-diamide-induced-TSS-Sbay-coords.tsv",
                tss_scer = "DSK2-diamide_TSS-seq-sense.tsv.gz",
                tss_smik = "smik-tss-seq.tsv.gz",
                tss_sbay = "sbay-tss-seq.tsv.gz",
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)

    df_peaks_scer = read_tsv(peak_scer) %>%
        filter(orf_name=="DSK2") %>%
        select(start, end, name, strand) %>%
        mutate(species="cerevisiae",
               annotation="peak")

    df_annotation = import_bed(transcripts_scer, "cerevisiae", "transcript") %>%
        bind_rows(import_bed(transcripts_smik, "mikatae", "transcript")) %>%
        bind_rows(import_bed(transcripts_sbay, "bayanus", "transcript")) %>%
        filter(name=="DSK2") %>%
        mutate(reference_point=start) %>%
        bind_rows(import_bed(orfs_scer, "cerevisiae", "orf")) %>%
        bind_rows(import_bed(orfs_smik, "mikatae", "orf")) %>%
        bind_rows(import_bed(orfs_sbay, "bayanus", "orf")) %>%
        filter(name=="DSK2") %>%
        bind_rows(import_bed(pace_scer, "cerevisiae", "pace")) %>%
        bind_rows(import_bed(pace_smik, "mikatae", "pace")) %>%
        bind_rows(import_bed(pace_sbay, "bayanus", "pace")) %>%
        bind_rows(read_tsv(peak_smik,
                           col_names = FALSE) %>%
                      filter(X16 %in% df_peaks_scer[["name"]]) %>%
                      select(start=X2,
                             end=X3,
                             name=X4,
                             strand=X6) %>%
                      mutate(species="mikatae",
                             annotation="peak")) %>%
        bind_rows(read_tsv(peak_sbay,
                           col_names = FALSE) %>%
                      filter(X16 %in% df_peaks_scer[["name"]]) %>%
                      select(start=X2,
                             end=X3,
                             name=X4,
                             strand=X6) %>%
                      mutate(species="bayanus",
                             annotation="peak")) %>%
        bind_rows(df_peaks_scer) %>%
        group_by(species) %>%
        mutate_at(vars(start, end),
                  ~((.-min(reference_point, na.rm=TRUE))/1000)) %>%
        select(-reference_point) %>%
        ungroup() %>%
        mutate(species = ordered(species,
                                 levels = c("cerevisiae", "mikatae", "bayanus"),
                                 labels = c("S. cerevisiae",
                                            "S. mikatae",
                                            "S. uvarum")),
               group = ordered("oxidative stress", levels=c("YPD", "oxidative stress"))) %>%
        left_join(., filter(., annotation=="orf") %>%
                      select(orf_end=end, species),
                  by="species")

    df_tss = read_tsv(tss_scer,
             col_names = c("group",
                           "sample",
                           "annotation",
                           "assay",
                           "index",
                           "position",
                           "signal")) %>%
        select(-c(annotation, index, assay)) %>%
        mutate(species="cerevisiae") %>%
        bind_rows(read_tsv(tss_smik,
                           col_names = c("group",
                                         "sample",
                                         "annotation",
                                         "index",
                                         "position",
                                         "signal")) %>%
                      select(-c(annotation, index)) %>%
                      mutate(species="mikatae")) %>%
        bind_rows(read_tsv(tss_sbay,
                           col_names = c("group",
                                         "sample",
                                         "annotation",
                                         "index",
                                         "position",
                                         "signal")) %>%
                      select(-c(annotation, index)) %>%
                      mutate(species="bayanus")) %>%
        group_by(group, position, species) %>%
        summarise(mean = mean(signal)) %>%
        ungroup() %>%
        mutate(group = case_when(group=="diamide" ~ "oxidative stress",
                                 group=="Sbay-diamide" ~ "oxidative stress",
                                 group=="Smik-diamide" ~ "oxidative stress",
                                 group=="untreated" ~ "YPD",
                                 group=="Sbay-YPD" ~ "YPD",
                                 group=="Smik-YPD" ~ "YPD"),
               group = ordered(group, levels = c("YPD", "oxidative stress")),
               species = ordered(species,
                                 levels = c("cerevisiae", "mikatae", "bayanus"),
                                 labels = c("S. cerevisiae",
                                            "S. mikatae",
                                            "S. uvarum"))) %>%
        group_by(species) %>%
        mutate(mean = scales::rescale(mean))

    figure_6a = ggplot() +
        geom_rect(data = df_annotation %>%
                      filter(annotation=="peak") %>%
                      select(-c(name, strand, annotation, orf_end)) %>%
                      complete(group, nesting(start, end, species)),
                  aes(xmin=start, xmax=end),
                  ymin=0, ymax=1.24,
                  fill="grey95") +
        geom_col(data = df_tss,
                 aes(x=position,
                     y=mean,
                     color=group),
                 size=0.2) +
        geom_segment(data = df_annotation %>%
                         filter(annotation=="transcript"),
                     aes(xend=end),
                     x=0,
                     y=1.07,
                     yend=1.07,
                     color="grey50") +
        geom_polygon(data = df_annotation %>%
                         filter(annotation=="orf") %>%
                         select(-c(name, strand, annotation, orf_end)) %>%
                         mutate(notch = 0.93*(end-start)+start) %>%
                         gather(key=point, value=position, -c(species, group)) %>%
                         mutate(y=case_when(point %in% c("start", "notch") ~ list(c(0.92, 1.22)),
                                            point=="end" ~ list(1.07))) %>%
                         unnest(y) %>%
                         group_by(species) %>%
                         arrange(position, y, .by_group=TRUE) %>%
                         mutate(point_order = c(1,2,5,3,4)) %>%
                         arrange(point_order, .by_group=TRUE),
                     aes(x=position, y=y),
                     fill="grey80") +
        geom_text(data = df_annotation %>%
                      filter(annotation=="orf"),
                  aes(x=(start+end*0.93)/2+start),
                  y=1.07,
                  size=10/72*25.4,
                  label="DSK2",
                  family="FreeSans",
                  fontface="italic") +
        geom_rect(data = df_annotation %>%
                      filter(annotation=="pace" &
                             end <= orf_end*0.93),
                  aes(xmin=start, xmax=end, fill=strand),
                  ymin=0.92, ymax=1.22) +
        geom_polygon(data = df_annotation %>%
                         filter(annotation=="pace" &
                                start > 0.93*orf_end &
                                end < orf_end) %>%
                     select(start, end, strand, orf_end, group, species) %>%
                     gather(key, x, -c(strand, orf_end, group, species)) %>%
                     mutate(y = -0.15/(0.07*orf_end)*x + 0.15/0.07,
                            yn = -y,
                            y = y+1.07,
                            yn = yn+1.07) %>%
                     gather(key, y, -c(strand, orf_end, key, x, group, species)) %>%
                     arrange(x, y) %>%
                     mutate(point_order = c(1,2,4,3)) %>%
                     arrange(point_order),
                     aes(x=x, y=y, fill=strand)) +
        geom_text(data = df_annotation %>%
                      distinct(species, group) %>%
                      # complete(species, group),
                      complete(species, group) %>%
                      filter(group=="YPD"),
                  aes(label=species),
                  x=-0.34, y=1,
                  hjust=0, vjust=1,
                  size=10/72*25.4,
                  family="FreeSans",
                  fontface="italic") +
        geom_text(data = df_annotation %>%
                      filter(species=="S. cerevisiae") %>%
                      distinct(species, group) %>%
                      complete(group, nesting(species)),
                  aes(label=group),
                  x=-0.34, y=0.65,
                  hjust=0, vjust=1,
                  size=10/72*25.4,
                  family="FreeSans") +
        facet_grid(species + group ~ .) +
        scale_x_continuous(expand = c(0,0),
                           breaks = scales::pretty_breaks(n=3),
                           name=NULL,
                           labels = function(x) (case_when(x==0 ~ "genic TSS",
                                                           near(x, 1) ~ "1 kb",
                                                           TRUE ~ as.character(x)))) +
        scale_y_continuous(limits = c(0, 1.24),
                           breaks = c(0, 1),
                           expand = c(0,0),
                           name="relative sense TSS-seq signal") +
        scale_color_viridis(discrete=TRUE, end=0.3, option="inferno") +
        scale_fill_manual(values = c("grey30", "grey60"),
                          guide=FALSE) +
        theme_default +
        theme(legend.position = "none",
              panel.border = element_blank(),
              panel.grid = element_blank(),
              axis.line.y = element_line(size=0.25, color="grey70"),
              axis.text.x = element_text(size=10),
              axis.text.y = element_text(size=8),
              axis.title.y = element_text(size=8),
              panel.spacing.y = unit(0, "mm"),
              plot.margin = margin(l=2, unit="pt"))


    ggplot2::ggsave(pdf_out,
                    plot=figure_6a,
                    width=fig_width,
                    height=fig_height,
                    units="in",
                    device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     transcripts_scer = snakemake@input[["transcripts_scer"]],
     transcripts_smik = snakemake@input[["transcripts_smik"]],
     transcripts_sbay = snakemake@input[["transcripts_sbay"]],
     orfs_scer = snakemake@input[["orfs_scer"]],
     orfs_smik = snakemake@input[["orfs_smik"]],
     orfs_sbay = snakemake@input[["orfs_sbay"]],
     pace_scer = snakemake@input[["pace_scer"]],
     pace_smik = snakemake@input[["pace_smik"]],
     pace_sbay = snakemake@input[["pace_sbay"]],
     peak_scer = snakemake@input[["peak_scer"]],
     peak_smik = snakemake@input[["peak_smik"]],
     peak_sbay = snakemake@input[["peak_sbay"]],
     tss_scer = snakemake@input[["tss_scer"]],
     tss_smik = snakemake@input[["tss_smik"]],
     tss_sbay = snakemake@input[["tss_sbay"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

