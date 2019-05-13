
extract_motif_ic = function(motif){
    motif %>%
        magrittr::extract("motif") %>%
        t() %>%
        as_tibble() %>%
        rowid_to_column(var="position") %>%
        gather(key="base",
               value="IC",
               -position) %>%
        group_by(position) %>%
        arrange(IC) %>%
        mutate(base_low=lag(cumsum(IC), default=0),
               base_high=base_low+IC) %>%
        left_join(ggseqlogo:::get_font("helvetica_bold"), by=c("base"="letter")) %>%
        group_by(position, base) %>%
        mutate(x=x+position-1.5,
               y=scales::rescale(y,
                                 to=c(first(base_low), first(base_high)))) %>%
        return()
}

parse_motifs = function(meme_path,
                        direction="sense",
                        category_id){
    motifs = universalmotif::read_meme(meme_path) %>%
        universalmotif::convert_type("ICM",
                                     relative_entropy=TRUE)
    if (direction == "antisense"){
        motifs %<>%
            universalmotif::motif_rc()
    }

    if (length(motifs) == 1){
        motifs %<>%
            list()
    }

    df = tibble()
    for (i in 1:length(motifs)){
        df = motifs[[i]] %>%
            extract_motif_ic() %>%
            mutate(motif=i) %>%
            bind_rows(df, .)
    }

    df %>%
        mutate(direction=direction,
               category=category_id) %>%
        return()

}

motif_plot = function(df,
                      df_annotation,
                      title){
    ggplot() +
        geom_hline(yintercept = 0,
                   size=0.25,
                   color="grey70") +
        geom_polygon(data=df,
                     aes(x=x,
                         y=y,
                         fill=base,
                         group=interaction(position,
                                           base,
                                           direction)),
                     alpha=0.95) +
        geom_text(data=df_annotation,
                  aes(label=paste0("E-value=", eval,
                                   "\nn=", nsites)),
                  x=29.5,
                  y=2,
                  hjust=1,
                  vjust=1,
                  size=5/72*25.4,
                  family="FreeSans") +
        facet_grid(motif~direction) +
        scale_fill_manual(values = c('#109648', '#255C99', '#F7B32B', '#D62839', '#D62839'),
                          breaks = c('A','C','G','T','U'),
                          guide=FALSE) +
        scale_y_continuous(limits = c(NA, 2),
                           breaks = scales::pretty_breaks(n=2),
                           expand=c(0,0),
                           name = "bits") +
        scale_x_continuous(expand = c(0,0),
                           limits = c(NA, 29.5),
                           breaks = scales::pretty_breaks(n=3),
                           labels = function(x) if_else(x==20,
                                                        paste(x, "nt"),
                                                        as.character(x))) +
        ggtitle(title) +
        theme_default +
        theme(axis.title.y = element_text(angle=0,
                                          hjust=1,
                                          vjust=0.5),
              axis.text.y = element_text(size=5),
              axis.title.x = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(size=0.25,
                                       color="grey70"),
              strip.text = element_text(color="black",
                                        size=7),
              strip.text.y = element_blank(),
              plot.title = element_text(size=7))
}

main = function(theme_spec,
                fonts_path,
                antisense_path,
                fig_width, fig_height,
                pdf_out){
    source(theme_spec)
    library(universalmotif)

    ttf_import(fonts_path)
    loadfonts()

    df = parse_motifs(antisense_path,
                 category_id="antisense") %>%
        bind_rows(parse_motifs(antisense_path,
                               direction="antisense",
                               category="antisense")) %>%
        mutate(category=ordered(category,
                                levels=c("intragenic",
                                         "antisense")),
               direction=ordered(direction,
                                 levels=c("sense",
                                          "antisense"),
                                 labels=c("motif",
                                          "reverse complement")))

get_motif_metadata = function(meme_path){
    motifs = universalmotif::read_meme(meme_path)

    if (length(motifs) == 1){
        motifs %<>%
            list()
    }

    df = tibble()
    for (i in 1:length(motifs)){
        df = motifs[[i]] %>%
            as.data.frame() %>%
            select(nsites,
                   eval) %>%
            mutate(motif=i) %>%
            bind_rows(df, .)
    }
    df %>%
        mutate(direction=ordered("motif",
                                 levels=c("motif",
                                          "reverse complement"))) %>%
        return()
}

    df_annotation = get_motif_metadata(antisense_path)

    motif_plots = motif_plot(df = df,
                             df_annotation= get_motif_metadata(antisense_path),
                             title="antisense")

    ggsave(pdf_out,
           plot=motif_plots,
           width=fig_width,
           height=fig_height,
           units="in",
           device=cairo_pdf)
}

main(theme_spec = snakemake@input[["theme"]],
     fonts_path = snakemake@input[["fonts_path"]],
     antisense_path = snakemake@input[["antisense"]],
     fig_width = snakemake@params[["width"]],
     fig_height = snakemake@params[["height"]],
     pdf_out = snakemake@output[["pdf"]])

