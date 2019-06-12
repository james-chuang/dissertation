
import = function(path, sample_list){
    read_tsv(path, col_names = c("group", "sample", "annotation", "assay", "index", "position", "signal")) %>%
        filter((sample %in% sample_list) & ! is.na(signal)) %>%
        mutate(signal=scales::rescale(signal),
               group = fct_inorder(group, ordered=TRUE)) %>%
        group_by(group, annotation, position) %>%
        summarise(mid = median(signal),
                  low = quantile(signal, 0.25),
                  high = quantile(signal, 0.75)) %>%
        ungroup() %>%
        return()
}

main = function(theme_spec = "thesis_theme.R",
                mnase_data = "transcripts-w-anti-plusonenucdistsort_MNase-seq.tsv.gz",
                tss_peaks = "depleted-v-non-depleted_tss-seq-spikenorm-peaks-diffexp-results-antisense-up.tsv",
                nuc_dyads = "S_pombe_alltranscripts-plus-one-nuc-dyad.bed",
                fig_width = 14,
                fig_height = 6.65,
                pdf_out = "test.pdf"){
    source(theme_spec)
    library(cowplot)

    sample_list = c("non-depleted-1, non-depleted-2", "non-depleted-3",
                    "depleted-1", "depleted-2", "depleted-3")

    df = import(mnase_data,
                sample_list=sample_list)

    bed = read_tsv(nuc_dyads,
                   col_names=c("chrom", "start", "end", "name", "score", "strand")) #%>%

    tss_df = read_tsv(tss_peaks) %>%
        left_join(bed, by=c("transcript_name"="name"),
                  suffix = c("_tss", "_nuc")) %>%
        group_by(name) %>%
        arrange(desc(score_nuc)) %>%
        slice(1) %>%
        ungroup() %>%
        mutate(nuc_to_tss_distance = if_else(transcript_strand=="+",
                                             start_tss + peak_summit - start_nuc,
                                             end_nuc - (start_tss + peak_summit + 1)))

    tss_df %>%
        arrange(nuc_to_tss_distance) %>%
        select(chrom_nuc, start_nuc, end_nuc, transcript_name, nuc_to_tss_distance, strand_nuc) %>%
        write_tsv("transcripts_w_antisense_sortbyplusonenucdist.bed", col_names=FALSE)

    mnase_plot = ggplot() +
        # geom_vline(xintercept=seq(0, 1.5, 0.152),
        #            color="grey90",
        #            linetype="dashed") +
        geom_ribbon(data = df,
                    aes(x=position, fill=group,
                        ymin=low, ymax=high),
                    alpha=0.2,
                    linetype='blank') +
        geom_line(data = df,
                    aes(x=position, color=group,
                        y=mid),
                  size=1,
                  alpha=0.85) +
        scale_x_continuous(breaks = seq(0, 1, 0.5),
                           labels = c("+1 dyad", "0.5", "1 kb"),
                           name = NULL,
                           expand = c(0,0)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n=2),
                           limits = c(0, max(df[["high"]])*1.05),
                           expand = c(0,0.002),
                           name = "normalized counts") +
        ggtitle("MNase-seq dyad signal") +
        scale_color_few() +
        scale_fill_few() +
        theme_default_presentation +
        theme(legend.key.height = unit(14, "pt"),
              legend.key.width = unit(18, "pt"),
              # axis.text.x = element_text(size=12,
              #                            margin=margin(t=1, unit="pt")),
              axis.text.x = element_blank(),
              axis.title.y = element_text(size=10),
              panel.grid = element_blank(),
              legend.justification=c(0.5, 1),
              legend.position = c(0.13, 0.96),
              legend.background = element_blank(),
              plot.margin=margin(b=4, unit="pt"))

    tss_plot = ggplot() +
        geom_vline(xintercept=seq(0, 1.5, 0.152),
                   color="grey50",
                   linetype="dashed") +
        geom_histogram(data= tss_df %>%
                           filter(nuc_to_tss_distance <= 1500),
                       aes(x=nuc_to_tss_distance / 1000,
                           y=..count..),
                 fill=viridis(1),
                 color=viridis(1),
                 size=0.1,
                 binwidth=0.01,
                 alpha=0.85) +
        scale_x_continuous(limits= c(-0.51,
                                     1.51),
                           breaks = seq(0, 1, 0.5),
                           labels = c("+1 dyad", "0.5", "1 kb"),
                           name = NULL,
                           expand = c(0,0)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n=3),
                           limits = c(0, 26),
                           expand = c(0,0.1),
                           name = "count") +
        ggtitle("antisense TSS frequency") +
        theme_default_presentation +
        theme(axis.text.x = element_text(size=12,
                                         margin=margin(t=1, unit="pt")),
              axis.title.y = element_text(size=10),
              plot.margin=margin(t=4, unit="pt"),
              panel.grid=element_blank())

    figure = plot_grid(mnase_plot,
                       tss_plot,
                       align="v",
                       axis="lr",
                       ncol=1,
                       rel_heights=c(0.55,0.45))

    ggsave(pdf_out,
           plot=figure,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     mnase_data = snakemake@input[["mnase_data"]],
     tss_peaks = snakemake@input[["tss_peaks"]],
     nuc_dyads = snakemake@input[["nuc_dyads"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

