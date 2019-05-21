import_bed = function(path){
    read_tsv(path,
             col_names = c("chrom", "start", "end",
                           "name", "score", "strand")) %>%
        select(-c(chrom, score)) %>%
        return()
}

main = function(theme_spec,
                mnase_chip_paths = list.files(".", "*reduced.tsv.gz"),
                tss_path = "DSK2_TSSseq.tsv.gz",
                tfiib_path ="DSK2_TFIIBnexus.tsv.gz",
                tfiib_diffbind_path = "diamide-v-YPD_tfiib-chipnexus-libsizenorm-peaks-diffbind-results-intragenic-all.tsv",
                transcript_annotation = "Scer_transcripts_w_verifiedORFs.bed",
                orf_annotation = "Scer_nondubious_ORFs_and_blocked_reading_frames.bed",
                dsk2_pace_annotation = "DSK2_PACE_core_sites.bed",
                gene_ids = "DSK2",
                up_distances = c(350),
                down_distances = c(0),
                modification ="h3k4me3",
                refpointlabel = "genic TSS",
                trim_pct,
                x_label_dist,
                fig_width, fig_height,
                svg_out, pdf_out, png_out, grob_out){
    source(theme_spec)
    library(cowplot)
    library(ggridges)

    dsk = read_tsv(transcript_annotation,
                   col_names = c("chrom", "start", "end", "name", "score", "strand")) %>%
        filter(name=="DSK2")

    pace_sites = read_tsv(dsk2_pace_annotation,
                      col_names = c("chrom", "start", "end",
                                    "name", "score", "strand")) %>%
        mutate_at(vars(start, end),
                  ~(.-dsk[["start"]]))

    tfiib_peak_annotation = read_tsv(tfiib_diffbind_path) %>%
        filter(orf_name=="DSK2") %>%
        mutate_at(vars(start, end),
                  ~(.-dsk[["start"]]))

    df = tibble()
    for (input in mnase_chip_paths){
        df = read_tsv(input) %>%
            bind_rows(df, .)
    }

    df %<>%
        select(-c(chip_low, chip_high, input_low, input_high,
                  enrichment_low, enrichment_high)) %>%
        drop_na() %>%
        mutate(log_enrichment = log2(enrichment_mid+1e-1),
               input_mid = scales::rescale(input_mid,
                                           from=c(0, df %>%
                                                      filter(! (group %in% c("h4k12ac", "h4k5ac", "h4k8ac"))) %>%
                                                      pull(input_mid) %>%
                                                      max(na.rm=TRUE))),
               input_mid = pmin(input_mid, 1),
               assay = "mnase_chip") %>%
        group_by(group) %>%
        mutate(log_enrichment = pmin(log_enrichment, max(log_enrichment[! is.infinite(log_enrichment)])),
               relative_log_enrichment =
                   (log_enrichment - mean(log_enrichment, na.rm=TRUE)) /
                   sd(log_enrichment, na.rm=TRUE))

    df_tfiib = read_tsv(tfiib_path,
                        col_names = c("assay", "time", "willzyx",
                                      "name", "index", "position", "input_mid")) %>%
        select(time, position, input_mid) %>%
        group_by(time, position) %>%
        summarise(input_mid = mean(input_mid)) %>%
        mutate(group = "",
               input_mid = scales::rescale(input_mid),
               assay = "tfiib_nexus") %>%
        complete(nesting(time, position, input_mid, assay), group = unique(df[["group"]])) %>%
        filter(group != "")
    df_tss = read_tsv(tss_path,
                        col_names = c("assay", "time", "willzyx",
                                      "name", "index", "position", "input_mid")) %>%
        select(time, position, input_mid) %>%
        group_by(time, position) %>%
        summarise(input_mid = mean(input_mid)) %>%
        mutate(group = "",
               input_mid = scales::rescale(input_mid),
               assay = "tss_seq") %>%
        complete(nesting(time, position, input_mid, assay), group = unique(df[["group"]])) %>%
        filter(group != "")
    df %<>%
        bind_rows(df_tfiib) %>%
        bind_rows(df_tss)

    ridgelines = ggplot() +
        annotate(geom="rect",
                 xmin=tfiib_peak_annotation[["start"]]/1000,
                 xmax=tfiib_peak_annotation[["end"]]/1000,
                 ymin=-10,
                 ymax=60,
                 fill="grey95") +
        geom_density_ridges_gradient(data = df %>%
                                         filter(group==modification),
                                     aes(x=position,
                                         y=time,
                                         height=input_mid,
                                         group=interaction(time, assay),
                                         fill=relative_log_enrichment,
                                         color=(assay=="tss_seq")),
                                     stat="identity",
                                     scale=1,
                                     panel_scaling=FALSE,
                                     size=0.2) +
        scale_fill_viridis(option="viridis",
                           name = "relative\nH3K4me3\nenrichment",
                           guide=guide_colorbar(barwidth=0.5,
                                                barheight=5,
                                                title.position="right",
                                                title.hjust=0,
                                                label.position="left"),
                           limits = c(quantile(df[["relative_log_enrichment"]], trim_pct, na.rm=TRUE),
                                      quantile(df %>%
                                                   filter(! is.infinite(relative_log_enrichment)) %>%
                                                   pull(relative_log_enrichment),
                                               1-trim_pct, na.rm=TRUE)),
                           oob = scales::squish,
                           breaks = scales::pretty_breaks(n=3),
                           na.value="#0000004c") +
        scale_color_manual(values = c("black", viridis(1)),
                           guide=FALSE) +
        scale_y_reverse(breaks = c(0, 4, 8, 15, 30, 45, 60),
                        sec.axis = sec_axis(~.,
                                            name = "smoothed MNase dyads,\nTFIIB protection,\nand sense TSS-seq",
                                            breaks = NULL),
                        expand = c(0,0),
                        name = "minutes in\noxidative\nstress") +
        scale_x_continuous(expand = c(0,0),
                           breaks = scales::pretty_breaks(n=3),
                           labels = function(x) case_when(x==0 ~ refpointlabel,
                                                          x==x_label_dist ~ paste(x, "kb"),
                                                          TRUE ~ as.character(x)),
                           name = NULL) +
        theme_default +
        theme(text = element_text(color="black", size=10),
              strip.background = element_blank(),
              strip.text = element_text(size=10, hjust=0),
              axis.ticks.y = element_blank(),
              axis.text.y = element_text(vjust=0, size=10),
              axis.text.x = element_text(size=10),
              axis.title.y.left = element_text(angle=0, hjust=1, vjust=0.5),
              axis.title.y.right = element_text(angle=0, hjust=0, vjust=0.8),
              panel.border = element_blank(),
              panel.grid.minor.y = element_blank(),
              legend.direction="vertical",
              legend.title = element_text(size=10),
              legend.text = element_text(size=8),
              legend.position=c(1.02, 0.299),
              legend.justification = c(0, 0.5))

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
        mutate(data_start = -up_distances,
               data_end = transcript_end + down_distances)

    gene_diagram = ggplot(data = annotation) +
        geom_segment(aes(x=data_start, xend=data_end),
                     y=0, yend=0,
                     size=0, color=NA) +
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
        geom_rect(data = pace_sites %>%
                      filter(end <= annotation[["orf_end"]]*0.93),
                  aes(xmin=start, xmax=end, fill=strand),
                  ymin=-1, ymax=1) +
        geom_polygon(data = pace_sites %>%
                         filter(end %>%
                                    between(annotation[["orf_end"]]*0.93,
                                            annotation[["orf_end"]])) %>%
                         select(start, end, strand) %>%
                         gather(key, x, -strand) %>%
                         mutate(y = -1/(0.07*annotation[["orf_end"]])*x + 1/0.07,
                                yn = -y) %>%
                         gather(key, y, -c(strand, key, x)) %>%
                         mutate(point_order = c(1,2,4,3)) %>%
                         arrange(point_order),
                     aes(x=x, y=y, fill=strand)) +
        geom_text(aes(x=(orf_start+orf_end)/2,
                      # label=paste0("italic(\"", gene_id, "\")")),
                      label=gene_id),
                  y=0, size=10/72*25.4, #parse=TRUE,
                  family="FreeSans",
                  fontface="italic") +
        scale_fill_manual(values=c("grey30", "grey60")) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(limits = c(-1.5, 1.5), expand=c(0,0)) +
        theme_void() +
        theme(strip.text.x = element_blank(),
              plot.margin = margin(b=-2, unit="pt"),
              legend.position="none")

    figure_4a = plot_grid(gene_diagram, ridgelines,
              align = "v", axis="lr",
              ncol=1,
              rel_heights = c(0.08, 1))

    ggplot2::ggsave(pdf_out,
                    plot=figure_4a,
                    width=fig_width,
                    height=fig_height,
                    units="in",
                    device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     mnase_chip_paths = snakemake@input[["mnase_chip_paths"]],
     tss_path = snakemake@input[["tss_path"]],
     tfiib_path = snakemake@input[["tfiib_path"]],
     tfiib_diffbind_path = snakemake@input[["tfiib_diffbind_path"]],
     transcript_annotation = snakemake@input[["transcript_annotation"]],
     orf_annotation = snakemake@input[["orf_annotation"]],
     dsk2_pace_annotation = snakemake@input[["dsk2_pace_annotation"]],
     gene_ids = snakemake@params[["gene_ids"]],
     up_distances = snakemake@params[["up_distances"]],
     down_distances = snakemake@params[["down_distances"]],
     modification = snakemake@params[["modification"]],
     refpointlabel = snakemake@params[["refpointlabel"]],
     trim_pct = snakemake@params[["trim_pct"]],
     x_label_dist = snakemake@params[["x_label_dist"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     svg_out = snakemake@output[["svg"]],
     pdf_out = snakemake@output[["pdf"]],
     png_out = snakemake@output[["png"]],
     grob_out = snakemake@output[["grob"]])

