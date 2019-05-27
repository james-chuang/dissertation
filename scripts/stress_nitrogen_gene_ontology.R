main = function(theme_spec = "thesis_theme.R",
                nitrogen_ontology_path = "nitrogen-v-SC_tfiib-chipnexus-libsizenorm-genic-up-gene-ontology-results.tsv",
                fig_width,
                fig_height,
                pdf_out){
    source(theme_spec)

    nitrogen_ontology = read_tsv(nitrogen_ontology_path) %>%
        mutate_at(vars(over_represented_pvalue,
                       under_represented_pvalue),
                  ~-log10(.)) %>%
        #mutate(term = case_when(#term=="cellular amino acid metabolic process" ~
        #                        #    "cellular amino acid\nmetabolic process",
        #                        #term=="transmembrane transporter activity" ~
        #                        #    "transmembrane\ntransporter activity",
        #                        term=="nucleobase-containing small molecule metabolic process" ~
        #                            "nucleobase-containing\nsmall molecule metabolic process",
        #                        TRUE ~ term)) %>%
        mutate(term = fct_inorder(term))

    ontology_fdr_cutoff = 0.01

    nitrogen_ontology_plot = ggplot(data = nitrogen_ontology %>%
                                        filter(over_represented_pvalue > -log10(ontology_fdr_cutoff))) +
        geom_hline(yintercept=-log10(ontology_fdr_cutoff),
                   linetype="dotted",
                   alpha=0.8) +
        geom_col(aes(x=fct_rev(term),
                     y=over_represented_pvalue,
                     fill=ontology),
                 alpha=0.9) +
        coord_flip() +
        scale_y_continuous(limits = c(0,
                                      1.05 * max(nitrogen_ontology[["over_represented_pvalue"]], na.rm=TRUE)),
                           breaks = scales::pretty_breaks(n=3),
                           expand = c(0,0),
                           name = expression("-log"[10] * "FDR")) +
        scale_fill_manual(values = ptol_pal()(4),
                          labels = c("biological process",
                                     "cellular compartment",
                                     "molecular function")) +
        ggtitle("enriched gene ontology terms",
                "nitrogen stress") +
        theme_default +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_text(size=10),
              legend.position = c(0.95, 0.45),
              legend.text = element_text(size=8),
              axis.ticks.y = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              plot.margin = margin(r=3, unit="pt"))

    ggplot2::ggsave(pdf_out,
                    plot=nitrogen_ontology_plot,
                    width=fig_width,
                    height=fig_height,
                    units="in",
                    device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     nitrogen_ontology_path = snakemake@input[["nitrogen_ontology"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

