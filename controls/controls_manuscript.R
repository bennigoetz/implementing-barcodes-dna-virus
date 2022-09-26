library(kableExtra)
library(magrittr)
## library(colorspace)

library(showtext)
library(patchwork)

library(ggpubr)

library(tidyverse)

font_add('LiberationMono',
         regular = '/usr/share/fonts/truetype/liberation2/LiberationMono-Regular.ttf',
         bold = '/usr/share/fonts/truetype/liberation2/LiberationMono-Bold.ttf')

font_add_google('IBM Plex Mono', 'IBM_Plex_Mono')

showtext_auto()

ten_control_barcodes = c(BC7_1='TCACAGGGGTAA', BC7_2='ACAAGACCGGAA',
                         BC7_3='ATATAGAGCTGT', BC7_4='ACATACCTGCTA',
                         BC7_5='GTGTCAGGCACA', BC7_6='TGCCACTCTAGC',
                         BC7_8='CTCGATTCACTC', BC7_9='GAACCCGTGGAA',
                         BC7_12='CTGTATATTTTA', BC7_15='GAAACCATGACA')

colors_ten_control_barcodes = colorspace::sequential_hcl(10, palette = 'Hawaii')
names(colors_ten_control_barcodes) = as.vector(ten_control_barcodes)

## with_bg_label_mapping = c(names(BC1_barcode), 'Other')
## names(with_bg_label_mapping) = names(with_bg_color_mapping)

## Reasonable to ask how well-distanced control barcodes are
adist(ten_control_barcodes) %>% as.vector() %>% unique() %>% sort()
## [1]  0  4  5  6  7  8  9 10

filenames_10diff3_clusters = c('10diff3_L0_clusters_dummy.tsv', list.files(pattern='^10diff3_mp_L[0-3]_clusters.tsv'))

counts_by_L_dist = filenames_10diff3_clusters %>%
    map(~ read_tsv(., col_names = c('barcode', str_extract(., 'L\\d'), 'elements'), col_types = 'cic')) %>%
    map(~ select(., -elements)) %>%
    reduce(full_join, by = 'barcode') %>%
    replace(is.na(.), 0) %>%
    pivot_longer(-barcode, names_to='L_dist', values_to='count') %>%
    group_by(L_dist) %>%
    mutate(max_count = max(count), frac = count/max_count)


# --- BC1 1-Plasmid -------------------------------------------------------

cluster_distances = c('L1', 'L2', 'L3')

one_plasmid_counts = expand.grid(library = c('BC1_102', 'BC1_103', 'BC1_104', 'BC1_105'),
                                 cluster = c('Raw', 'L1', 'L2', 'L3'),
                                 stringsAsFactors = FALSE
                                 ) %>%
    ## filter(cluster != 'raw') %>%
    mutate(data = map2(.x = library, .y = cluster,
                       ~ if (.y != 'Raw') {
                             read_tsv(paste(.x, 'mp', .y, 'clusters.tsv', sep = '_'),
                                      col_names = c('barcode', 'count', 'elements'),
                                      progress = FALSE) %>%
                                 select(-elements)
                         }
                         else {
                             read_table(paste0('../barcodes_e0.1/JA19375/', .x, '_cutadapt_counts.txt'),
                                        col_names = c('count', 'barcode'),
                                        progress = FALSE)
                               })) %>%
    unnest(., data)

## Calculate what percentage of barcodes aren't the canonical barcode
one_plasmid_counts %>%
    group_by(library, cluster) %>%
    summarize(non_canonical_pct =
                  sum((barcode != ten_control_barcodes[['BC7_1']]) * count) /
                  sum(count)
              ) %>%
    mutate(non_canonical_pct = scales::percent(non_canonical_pct, accuracy = 0.01)) %>%
    arrange(match(cluster, c('Raw', 'L1', 'L2', 'L3')), .by_group = TRUE) %>%
    kbl(caption = 'Percentage Of Non-Canonical Barcodes') %>%
    column_spec(1, bold = TRUE) %>%
    collapse_rows(columns = 1, valign = 'top') %>%
    kable_styling(bootstrap_options='striped', full_width=FALSE) %>%
    save_kable('../plots/controls/tables/one_plasmid_pct_non_canonical_bcs.html')


one_plasmid_counts %>%
    filter(cluster == 'Raw') %>%
    group_by(library) %>%
    mutate(frac = count / sum(count)) %>%
    arrange(library, desc(frac)) %>%
    ungroup() %>%
    group_by(barcode) %>%
    mutate(max_frac_across_libraries = max(frac)) %>%
    ungroup() %>%
    filter(barcode %in% (top_n(., 60, wt = max_frac_across_libraries) %>%
                         pull(barcode) %>%
                         unique())
           ) %>%
    arrange(desc(max_frac_across_libraries)) %>%
    mutate(barcode = fct_reorder(as.factor(barcode), max_frac_across_libraries)) %>%
    left_join(tibble(library = c('BC1_102', 'BC1_103', 'BC1_104', 'BC1_105'),
                     library_name = c('1P-A', '1P-B', '1P-C', '1P-D'))  # This is to rename libraries
        ) %>%
    ggplot(aes(x = barcode, y = frac)) +
    geom_col(aes(fill = barcode)) +
    coord_flip() +
    facet_grid(cols = vars(library_name)) +
    ## labs(title = 'One-Plasmid Controls Unclustered', 
    ##      subtitle = 'Fraction Of Total Reads',
    ##      x = 'Barcode') +
    ## scale_fill_manual(name = 'Control',
    ##                   values = colors_ten_control_barcodes,
    ##                   limits = c(ten_control_barcodes[['BC7_1']], 'Other'),
    ##                   na.value = 'gray85') +
    labs(x = 'Barcodes') +
    scale_fill_manual(name = 'Barcode',
                      values = c(colors_ten_control_barcodes[[1]], 'gray85'),
                      limits = c(ten_control_barcodes[['BC7_1']], 'Error'),
                      na.value = 'gray85') +
    scale_y_continuous(labels = scales::label_percent(),
                       minor_breaks = NULL) +
    theme_minimal() +
    theme(panel.spacing.y = unit(2, "lines"),
          ## axis.title.x = element_blank(),
          axis.text.y = element_text(family = 'IBM_Plex_Mono', size = 10),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          axis.title.x = element_text(vjust = -1, size = 12),
          strip.text = element_text(size = 12),
          legend.text = element_text(family = 'IBM_Plex_Mono'),
          ## legend.position = 'none',
          plot.caption = element_text(hjust = 0, size = 14)
          ) -> p_one_plasmid_pct_total_bar

ggsave('one_plasmid_pct_unclustered_bar.pdf',
       path = '../plots/controls',
       plot = p_one_plasmid_pct_total_bar,
       width = 14)
       
ggsave('one_plasmid_pct_unclustered_sqrt_bar.pdf',
       plot = p_one_plasmid_pct_total_bar +
           ## scale_y_continuous('Percent Of Total Barcodes',
           ##                    trans = 'sqrt',
           ##                    labels = scales::label_percent(),
           ##                    minor_breaks = NULL) +
           scale_y_continuous('Square Root-Transformed Percentage of Total Read Counts',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              n.breaks = 3,
                              minor_breaks = NULL),
       path='../plots/controls',
       width=14)

# EPS Version
ggsave('one_plasmid_pct_unclustered_sqrt_bar.eps',
       plot = p_one_plasmid_pct_total_bar +
           ## scale_y_continuous('Percent Of Total Barcodes',
           ##                    trans = 'sqrt',
           ##                    labels = scales::label_percent(),
           ##                    minor_breaks = NULL) +
           scale_y_continuous('Square Root-Transformed Percentage Of Total Read Counts',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              minor_breaks = NULL),
       path='../plots/controls',
       device = 'eps',
       width=14)

# TIFF Version
ggsave('one_plasmid_pct_unclustered_sqrt_bar.tiff',
       plot = p_one_plasmid_pct_total_bar +
           ## scale_y_continuous('Percent Of Total Barcodes',
           ##                    trans = 'sqrt',
           ##                    labels = scales::label_percent(),
           ##                    minor_breaks = NULL) +
           scale_y_continuous('Square Root-Transformed Percentage Of Total Read Counts',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              minor_breaks = NULL),
       path='../plots/controls',
       device = 'tiff',
       dpi = 300,
       compression = 'lzw',
       width = 2200,
       units = 'px')

## -- BC1_103 -------------------------------------------------------------

one_plasmid_counts %>%
    filter(library == 'BC1_103') %>%
    group_by(cluster) %>%
    mutate(frac = count / sum(count)) %>%
    arrange(cluster, desc(frac)) %>%
    ungroup() %>%
    group_by(barcode) %>%
    mutate(max_frac_across_clusters = max(frac)) %>%
    ungroup() %>%
    filter(barcode %in% (group_by(., barcode) %>%
                         filter(frac == max(frac)) %>%
                         ungroup() %>%
                         top_n(15, wt = frac) %>%
                         pull(barcode))
           ) %>%
    complete(library, cluster, nesting(barcode), fill = list(count = 0, frac = 0, max_frac_across_clusters = 0)) %>%
    ## arrange(desc(max_frac_across_clusters)) %>%
    ## mutate(barcode = fct_inorder(as.factor(barcode), ordered = FALSE)) %>%
    ggplot(aes(x = reorder(barcode, max_frac_across_clusters), y = frac)) +
    geom_col(aes(fill = barcode)) +
    coord_flip() +
    facet_grid(cols = vars(factor(cluster, levels = c('Raw', 'L1', 'L2', 'L3')))) +
    labs(title = 'One-Plasmid Control BC1_103', 
         subtitle = 'Fraction Of Total Reads',
         x = 'Barcode') +
    scale_fill_manual(name = 'Control',
                      values = colors_ten_control_barcodes,
                      limits = c(ten_control_barcodes[['BC7_1']], 'Other'),
                      na.value = 'gray85') +
    scale_y_continuous('Percent Of Total Barcodes',
                       labels = scales::label_percent(),
                       minor_breaks = NULL) +
    theme_minimal() +
    theme(panel.spacing.y = unit(2, "lines"),
          axis.text.y = element_text(family = 'IBM_Plex_Mono'),
          axis.title = element_blank(),
          legend.text = element_text(family = 'IBM_Plex_Mono'),
          text = element_text(family = 'Arial')
          ) -> p_one_plasmid_BC1_103_Lparams_pct_total_bar

ggsave('one_plasmid_BC1_103_pct_clustered_bar.pdf',
       path = '../plots/controls',
       plot = p_one_plasmid_BC1_103_Lparams_pct_total_bar,
       width = 14)
       
ggsave('one_plasmid_BC1_103_pct_clustered_sqrt_bar.pdf',
       plot = last_plot() +
           scale_y_continuous('Percent Of Total Barcodes',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              minor_breaks = NULL) +
           labs(subtitle='Fraction Of Total Reads, Square Root'),
       path='../plots/controls',
       width=14)

# EPS/Arial Version
ggsave('one_plasmid_BC1_103_pct_clustered_sqrt_bar.eps',
       plot = last_plot() +
           scale_y_continuous('Percent Of Total Barcodes',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              minor_breaks = NULL) +
           labs(subtitle='Fraction Of Total Reads, Square Root'),
       path='../plots/controls',
       device = 'eps',
       width=14)

# --- 10diff1 Controls ----------------------------------------------------

diff_plasmid_counts = expand.grid(library = c('10diff1', '10diff2', '10diff3'),
                                  cluster = c('Raw', 'L1', 'L2', 'L3'),
                                  stringsAsFactors = FALSE) %>%
    ## filter(cluster != 'Raw') %>%
    mutate(data = map2(.x = library, .y = cluster,
                       ~ if (.y != 'Raw') {
                             read_tsv(paste(.x, 'mp', .y, 'clusters.tsv', sep = '_'),
                                      col_names = c('barcode', 'count', 'elements'),
                                      progress = FALSE) %>%
                                 select(-elements)
                         }
                         else {
                             read.table(paste0('../barcodes_e0.1/JA19375/', .x, '_cutadapt_counts.txt'),
                                              header = FALSE,
                                              col.names = c('count', 'barcode'),
                                              stringsAsFactors = FALSE)
                               })) %>%
    unnest(., data)

### 10diff1 (10P-A)

diff_plasmid_counts %>%
    filter(library == '10diff3') %>%
    group_by(cluster) %>%
    mutate(frac = count / sum(count)) %>%
    arrange(cluster, desc(frac)) %>%
    ungroup() %>%
    group_by(barcode) %>%
    mutate(max_frac_across_clusters = max(frac)) %>%
    ungroup() %>%
    ## filter(barcode %in% (group_by(., barcode) %>%
    ##                      filter(frac == max(frac)) %>%
    ##                      ungroup() %>%
    ##                      top_n(15, wt = frac) %>%
    ##                      pull(barcode))
    ##        ) %>%
    filter(barcode %in% union(as.vector(ten_control_barcodes),
                                        group_by(., barcode) %>%
                                        filter(frac == max(frac)) %>%
                                        ungroup() %>%
                                        filter(!barcode %in% ten_control_barcodes) %>%
                                        top_n(5, wt = frac) %>%
                                        pull(barcode)
                              )
           ) %>%
    complete(library, cluster, nesting(barcode), fill = list(count = 0, frac = 0, max_frac_across_clusters = 0)) %>%
    ggplot(aes(x = reorder(barcode, max_frac_across_clusters), y = frac)) +
    geom_col(aes(fill = barcode)) +
    coord_flip() +
    facet_grid(cols = vars(factor(cluster, levels = c('Raw', 'L1', 'L2', 'L3')))) +
    labs(x = 'Barcode') +
    scale_fill_manual(name = 'Controls',
                      values = c(colors_ten_control_barcodes, 'gray85'),
                      limits = c(ten_control_barcodes, 'Other'),
                      na.value = 'gray85') +
    scale_y_continuous('Percent Of Total Barcodes',
                       labels = scales::label_percent(accuracy = 0.1),
                       breaks = c(0.02, 0.1, 0.25, 0.5),
                       minor_breaks = NULL) +
    theme_minimal() +
    theme(panel.spacing.y = unit(2, "lines"),
          axis.title = element_blank(),
          axis.text.y = element_text(family = 'IBM_Plex_Mono'),
          legend.text = element_text(family = 'IBM_Plex_Mono'),
          legend.position = 'none'
          ) -> p_diff_plasmid_10diff3_Lparams_pct_total_bar

ggsave('diff_plasmid_10diff3_pct_clustered_sqrt_bar.pdf',
       plot = p_diff_plasmid_10diff3_Lparams_pct_total_bar +
           scale_y_continuous('Percent Of Total Barcodes',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              breaks = c(.125, 0.25, 0.375, 0.5),
                              minor_breaks = NULL),
       path='../plots/controls',
       width=14)


diff_plasmid_counts %>%
    filter(library == '10diff2') %>%
    group_by(cluster) %>%
    mutate(frac = count / sum(count)) %>%
    arrange(cluster, desc(frac)) %>%
    ungroup() %>%
    group_by(barcode) %>%
    mutate(max_frac_across_clusters = max(frac)) %>%
    ungroup() %>%
    ## filter(barcode %in% (group_by(., barcode) %>%
    ##                      filter(frac == max(frac)) %>%
    ##                      ungroup() %>%
    ##                      top_n(15, wt = frac) %>%
    ##                      pull(barcode))
    ##        ) %>%
    filter(barcode %in% union(as.vector(ten_control_barcodes),
                                        group_by(., barcode) %>%
                                        filter(frac == max(frac)) %>%
                                        ungroup() %>%
                                        filter(!barcode %in% ten_control_barcodes) %>%
                                        top_n(5, wt = frac) %>%
                                        pull(barcode)
                              )
           ) %>%
    complete(library, cluster, nesting(barcode), fill = list(count = 0, frac = 0, max_frac_across_clusters = 0)) %>%
    ggplot(aes(x = reorder(barcode, max_frac_across_clusters), y = frac)) +
    geom_col(aes(fill = barcode)) +
    coord_flip() +
    facet_grid(cols = vars(factor(cluster, levels = c('Raw', 'L1', 'L2', 'L3')))) +
    labs(x = 'Barcode') +
    scale_fill_manual(name = 'Controls',
                      values = c(colors_ten_control_barcodes, 'gray85'),
                      limits = c(ten_control_barcodes, 'Other'),
                      na.value = 'gray85') +
    scale_y_continuous('Percent Of Total Barcodes',
                       labels = scales::label_percent(accuracy = 0.1),
                       breaks = c(0.02, 0.1, 0.25, 0.5),
                       minor_breaks = NULL) +
    theme_minimal() +
    theme(panel.spacing.y = unit(2, "lines"),
          axis.title = element_blank(),
          axis.text.y = element_text(family = 'IBM_Plex_Mono'),
          legend.text = element_text(family = 'IBM_Plex_Mono'),
          legend.position = 'none'
          ) -> p_diff_plasmid_10diff2_Lparams_pct_total_bar

ggsave('diff_plasmid_10diff2_pct_clustered_bar.pdf',
       path = '../plots/controls',
       plot = p_diff_plasmid_10diff2_Lparams_pct_total_bar,
       width = 14)
       
ggsave('diff_plasmid_10diff2_pct_clustered_sqrt_bar.pdf',
       plot = p_diff_plasmid_10diff2_Lparams_pct_total_bar +
           scale_y_continuous('Percent Of Total Barcodes',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              breaks = c(.125, 0.25, 0.375, 0.5),
                              minor_breaks = NULL),
       path='../plots/controls',
       width=14)


# --- Ten Plasmid Controls ------------------------------------------------

ten_plasmid_counts = expand.grid(library = c('1_01', '1_12', '10_10', '10_12'),
                                 cluster = c('Raw', 'L1', 'L2', 'L3'),
                                 stringsAsFactors = FALSE) %>%
    ## filter(cluster != 'Raw') %>%
    mutate(data = map2(.x = library, .y = cluster,
                       ~ if (.y != 'Raw') {
                             read_tsv(paste(.x, 'mp', .y, 'clusters.tsv', sep = '_'),
                                      col_names = c('barcode', 'count', 'elements'),
                                      progress = FALSE) %>%
                                 select(-elements)
                         }
                         else {
                             read.table(paste0('../barcodes_e0.1/JA19375/', .x, '_cutadapt_counts.txt'),
                                        header = FALSE,
                                        col.names = c('count', 'barcode'),
                                        fill = TRUE,  # for empty barcodes
                                        stringsAsFactors = FALSE)
                               })) %>%
    unnest(., data)

ten_plasmid_counts %>%
    filter(library == '1_12') %>%
    group_by(cluster) %>%
    mutate(frac = count / sum(count)) %>%
    arrange(cluster, desc(frac)) %>%
    ungroup() %>%
    group_by(barcode) %>%
    mutate(max_frac_across_clusters = max(frac)) %>%
    ungroup() %>%
    ## filter(barcode %in% (group_by(., barcode) %>%
    ##                      filter(frac == max(frac)) %>%
    ##                      ungroup() %>%
    ##                      top_n(15, wt = frac) %>%
    ##                      pull(barcode))
    ##        ) %>%
    filter(barcode %in% union(as.vector(ten_control_barcodes),
                                        group_by(., barcode) %>%
                                        filter(frac == max(frac)) %>%
                                        ungroup() %>%
                                        filter(!barcode %in% ten_control_barcodes) %>%
                                        top_n(5, wt = frac) %>%
                                        pull(barcode)
                              )
           ) %>%
    complete(library, cluster, nesting(barcode), fill = list(count = 0, frac = 0, max_frac_across_clusters = 0)) %>%
    ggplot(aes(x = reorder(barcode, max_frac_across_clusters), y = frac)) +
    geom_col(aes(fill = barcode)) +
    coord_flip() +
    facet_grid(cols = vars(factor(cluster, levels = c('Raw', 'L1', 'L2', 'L3')))) +
    ## labs(title = 'Ten??? Control 1_12',
    ##      subtitle = 'Fraction Of Total Reads',
    ##      x = 'Barcode') +
    scale_fill_manual(name = 'Controls',
                      values = c(colors_ten_control_barcodes, 'gray85'),
                      limits = c(ten_control_barcodes, 'Other'),
                      na.value = 'gray85') +
    scale_y_continuous('Percent Of Total Barcodes',
                       labels = scales::label_percent(accuracy = 1),
                       ## breaks = c(.125, 0.25, 0.375, 0.5),
                       minor_breaks = NULL) +
    theme_minimal() +
    theme(panel.spacing.y = unit(2, "lines"),
          axis.title = element_blank(),
          axis.text.y = element_text(family = 'IBM_Plex_Mono'),
          legend.text = element_text(family = 'IBM_Plex_Mono'),
          legend.position = 'none'
          ) -> p_ten_plasmid_1_12_Lparams_pct_total_bar

ggsave('ten_plasmid_1_12_pct_clustered_bar.pdf',
       path = '../plots/controls',
       plot = p_ten_plasmid_1_12_Lparams_pct_total_bar,
       width = 14)
       
ggsave('ten_plasmid_1_12_pct_clustered_sqrt_bar.pdf',
       plot = last_plot() +
           scale_y_continuous('Percent Of Total Barcodes',
                              trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              ## breaks = c(.125, 0.25, 0.375, 0.5),
                              minor_breaks = NULL),
       path='../plots/controls',
       width=14)


# --- Ten-Plasmid Controls All Together -----------------------------------

rename_10P = c(
    '10diff1' = '10P-A',
    '10diff2' = '10P-B',
    '10diff3' = '10P-C',
    '1_01' = '10P-D',
    '1_12' = '10P-E',
    '10_10' = '10P-F',
    '10_12' = '10P-G'
    )

all_ten_plasmids =
    bind_rows(diff_plasmid_counts, ten_plasmid_counts) %>%
    mutate(library = rename_10P[library])

plot_all_ten = all_ten_plasmids %>%
    group_by(library, cluster) %>%
    mutate(frac = count / sum(count)) %>%
    ## arrange(library, cluster, desc(frac)) %>%
    ## ungroup(cluster) %>%
    ## group_by(library, barcode) %>%
    ## ## mutate(max_frac_across_clusters = max(frac)) %>%
    ungroup() %>%
    ## group_by(library) %>%
    ## filter(barcode %in% ten_control_barcodes) %>%
    ## filter(barcode %in% (
    ##                                     group_by(., barcode) %>%
    ##                                     filter(frac == max(frac)) %>%
    ##                                     ungroup() %>%
    ##                                     filter(!barcode %in% ten_control_barcodes) %>%
    ##                                     top_n(5, wt = frac) %>%
    ##                                     pull(barcode)
    ##                           )
    ##        ) %>%
    filter(barcode %in% 
           (
               ## filter(., ! barcode %in% ten_control_barcodes) %>%
     group_by(., library, barcode) %>%
    top_n(1, wt = frac) %>%
    group_by(library) %>%
    top_n(5, wt = frac) %>%
    pull(barcode)
    )
    )
    

    
    filter(barcode %in% union(as.vector(ten_control_barcodes),
                                        group_by(., barcode) %>%
                                        filter(frac == max(frac)) %>%
                                        ungroup() %>%
                                        group_by(library) %>%
                                        filter(!barcode %in% ten_control_barcodes) %>%
                                        top_n(5, wt = frac) %>%
                                        pull(barcode) 
                              )
           )

## attempt with nesting from the beginning

## selected_barcodes = all_ten_plasmids %>%
##     group_by(library, cluster) %>%
##     mutate(frac = count / sum(count)) %>%
##     ungroup() %>%
##     split(.$library) %>%
##     imap(~ select(.x, - library)) %>%
##     imap(~ filter(.x, ! barcode %in% ten_control_barcodes)) %>%
##     imap(~ group_by(.x, barcode)) %>%
##     imap(~ mutate(.x, selection_metric = max(frac))) %>%
##     ## imap(~ mutate(.x, selection_metric = median(frac))) %>%
##     ## imap(~ mutate(.x, selection_metric = mean(frac))) %>%
##     imap(~ distinct(.x, barcode, .keep_all = TRUE)) %>%
##     imap(~ ungroup(.x)) %>%
##     imap(~ top_n(.x, 5, wt = selection_metric)) %>%
##     imap(~ pull(.x, barcode)) %>%
##     imap(~ union(.x, ten_control_barcodes)) %>%
##     enframe(name = 'library', value = 'selected_barcodes') %>%
##     unnest(cols = selected_barcodes)

selected_barcodes = all_ten_plasmids %>%
    group_by(library, cluster) %>%
    mutate(frac = count / sum(count)) %>%
    ungroup() %>%
    split(.$library) %>%
    imap(~ select(.x, - library)) %>%
    imap(~ filter(.x, ! barcode %in% ten_control_barcodes)) %>%
    imap(~ group_by(.x, cluster)) %>%
    imap(~ top_n(.x, 2, wt = frac)) %>%
    imap(~ pull(.x, barcode)) %>%
    imap(~ union(.x, ten_control_barcodes)) %>%
    enframe(name = 'library', value = 'selected_barcodes') %>%
    unnest(cols = selected_barcodes)

main_text_10P_controls = c('10P-A', '10P-D', '10P-G')

all_ten_plasmids %>%
    complete(library, cluster, barcode, fill = list(count = 0)) %>%
    group_by(library, cluster) %>%
    mutate(frac = count / sum(count)) %>%
    ungroup() %>%
    right_join(selected_barcodes %>% filter(library %in% main_text_10P_controls),
               by = c('library' = 'library', 'barcode' = 'selected_barcodes')) %>% 
    group_by(library) %>%
    nest() %>%
    mutate(plot = map(data, ~ ggplot(data = .x, aes(x = reorder(barcode, frac), y = frac)) +
                                geom_col(aes(fill = barcode)) +
                                coord_flip() +
                                facet_grid(cols = vars(factor(cluster, levels = c('Raw', 'L1', 'L2', 'L3')))) +
                                scale_fill_manual(name = 'Controls',
                                                  values = c(colors_ten_control_barcodes, 'gray85'),
                                                  limits = c(ten_control_barcodes, 'Other'),
                                                  na.value = 'gray85') +
                                scale_y_continuous('Percent Of Total Barcodes',
                                                   labels = scales::label_percent(accuracy = 1),
                                                   ## breaks = c(.125, 0.25, 0.375, 0.5),
                                                   minor_breaks = NULL) +
                                labs(tag = library) +
                                theme_minimal() +
                                theme(panel.spacing.y = unit(2, "lines"),
                                      axis.title = element_blank(),
                                      axis.text.y = element_text(family = 'IBM_Plex_Mono'),
                                      legend.text = element_text(family = 'IBM_Plex_Mono'),
                                      legend.position = 'none'
                                )
                      )
           ) %>%
    pull(plot) %>%
    reduce(`/`) +
    ## plot_annotation(
    ##     tag_levels = 'A',
    ##     tag_prefix = 'Fig. ',
    ##     theme = theme(plot.title = element_text(size = 18))
    ## ) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    ## xlab(label = 'Percent Of Total Barcodes') ->
    p_10P_ADG_L_params_pct_bar
    ## p_10P_A_D_G_L_params_pct_bar
    ## p_10P_C_D_L_params_pct_bar
    ## p_10P_L_params_pct_bar

ten_plasmid_L_plots = all_ten_plasmids %>%
    complete(library, cluster, barcode, fill = list(count = 0)) %>%
    group_by(library, cluster) %>%
    mutate(frac = count / sum(count)) %>%
    ungroup() %>%
    semi_join(selected_barcodes, by = c('library' = 'library', 'barcode' = 'selected_barcodes')) %>% 
    group_by(library) %>%
    nest() %>%
    mutate(plot = map(data, ~ ggplot(data = .x, aes(x = reorder(barcode, frac), y = frac)) +
                                geom_col(aes(fill = barcode)) +
                                coord_flip() +
                                facet_grid(cols = vars(factor(cluster, levels = c('Raw', 'L1', 'L2', 'L3')))) +
                                scale_fill_manual(name = 'Controls',
                                                  values = c(colors_ten_control_barcodes, 'gray85'),
                                                  limits = c(ten_control_barcodes, 'Other'),
                                                  na.value = 'gray85') +
                                scale_y_continuous('Percent of Total Barcodes',
                                                   labels = scales::label_percent(accuracy = 1),
                                                   ## breaks = c(.125, 0.25, 0.375, 0.5),
                                                   minor_breaks = NULL) +
                                labs(tag = library) +
                                theme_minimal() +
                                theme(panel.spacing.y = unit(2, "lines"),
                                      ## axis.title = element_blank(),
                                      ## axis.title.x = element_blank(),
                                      axis.title.y = element_text(vjust = 10),
                                      axis.text.y = element_text(family = 'IBM_Plex_Mono', size = 10),
                                      axis.title.x = element_text(vjust = -1),
                                      axis.text.x = element_text(size = 10),
                                      legend.text = element_text(family = 'IBM_Plex_Mono'),
                                      legend.position = 'none'
                                )
                      )
           ) %>%
    ungroup() %>%
    select(-data) %>%
    deframe()

p_10P_ADG_bar = ten_plasmid_L_plots[c('10P-A', '10P-D', '10P-G')]
p_10P_ADG_bar[['10P-A']] = p_10P_ADG_bar[['10P-A']] +
    theme(axis.title.x = element_blank()) #+
    ## theme(axis.title.y = element_blank())
p_10P_ADG_bar[['10P-D']] = p_10P_ADG_bar[['10P-D']] +
    theme(axis.title.x = element_blank())
## p_10P_ADG_bar[['10P-G']] = p_10P_ADG_bar[['10P-G']] +
##     theme(axis.title.y = element_blank())

p_10P_ADG_bar %<>%
    reduce(`/`) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12))

ggsave('controls_10P_ADG_bar.pdf',
       path = '../plots/controls',
       plot = p_10P_ADG_bar &
           scale_y_continuous('Square Root-Transformed Percentage of Total Read Counts',
                              trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL) &
           scale_x_discrete('Barcodes'),
       width = 15,
       height = 15
       )

p_10P_BCEF_bar = ten_plasmid_L_plots[c('10P-B', '10P-C', '10P-E', '10P-F')]
p_10P_BCEF_bar[['10P-B']] = p_10P_BCEF_bar[['10P-B']] +
    theme(axis.title.x = element_blank()) #+
    ## theme(axis.title.y = element_blank())
p_10P_BCEF_bar[['10P-C']] = p_10P_BCEF_bar[['10P-C']] +
    theme(axis.title.x = element_blank())
p_10P_BCEF_bar[['10P-E']] = p_10P_BCEF_bar[['10P-E']] +
    theme(axis.title.x = element_blank())
## p_10P_BCEF_bar[['10P-E']] = p_10P_BCEF_bar[['10P-E']] +
##     theme(axis.title.y = element_blank())

p_10P_BCEF_bar %<>%
    reduce(`/`) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12))

ggsave('controls_10P_BCEF_bar.pdf',
       path = '../plots/controls',
       plot = p_10P_BCEF_bar &
           scale_y_continuous('Square Root-Transformed Percentage of Total Read Counts',
                              trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL) &
           scale_x_discrete('Barcodes'),
       width = 15,
       height = 15
       )


ggsave('controls_all_10P_bar.pdf',
       path = '../plots/controls',
       plot = p_10P_L_params_pct_bar &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 30)

ggsave('controls_10P-C_10P-D_bar.pdf',
       path = '../plots/controls',
       plot = p_10P_C_D_L_params_pct_bar &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 10)

ggsave('controls_10P_ADG_bar.pdf',
       path = '../plots/controls',
       plot = p_10P_ADG_L_params_pct_bar +
           xlab(label = 'Percent Of Total Barcodes') &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 15)

ggsave('controls_10P-B-C-E-F_bar.pdf',
       path = '../plots/controls',
       plot = p_10P_B_C_E_F_L_params_pct_bar &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 20)



# --- Compose Single Sequence Of L Of Three Types Of Controls -------------

p_controls_BC1_103_10diff2_1_12_Lparams_pct_total_bar =
    (
        p_one_plasmid_BC1_103_Lparams_pct_total_bar +
        theme(legend.position = 'None')
    ) /
        p_diff_plasmid_10diff2_Lparams_pct_total_bar / 
        p_ten_plasmid_1_12_Lparams_pct_total_bar +
    plot_annotation(
        ## title = 'Percent Of Total Reads For Controls',
        ## subtitle = 'Square Root Transformed',
        tag_levels = 'A',
        tag_prefix = 'Fig. ',
        theme = theme(plot.title = element_text(size = 18))
        ) +
    plot_layout(guides = 'collect') &
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12))

p_controls_10diff2_1_12_Lparams_pct_total_bar =
        p_diff_plasmid_10diff2_Lparams_pct_total_bar / 
        p_ten_plasmid_1_12_Lparams_pct_total_bar +
    plot_annotation(
        ## title = 'Percent Of Total Reads For Controls',
        ## subtitle = 'Square Root Transformed',
        tag_levels = 'A',
        tag_prefix = 'Fig. ',
        theme = theme(plot.title = element_text(size = 18))
        ) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12))

p_controls_10diff3_1_12_Lparams_pct_total_bar =
        p_diff_plasmid_10diff3_Lparams_pct_total_bar / 
        p_ten_plasmid_1_12_Lparams_pct_total_bar +
    plot_annotation(
        ## title = 'Percent Of Total Reads For Controls',
        ## subtitle = 'Square Root Transformed',
        tag_levels = 'A',
        tag_prefix = 'Fig. ',
        theme = theme(plot.title = element_text(size = 18))
        ) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12))

ggsave('controls_BC1_103_10diff2_1_12_clustered_sqrt_bar.pdf',
       path = '../plots/controls',
       plot = p_controls_BC1_103_10diff2_1_12_Lparams_pct_total_bar &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 18)


ggsave('controls_10diff3_1_12_clustered_sqrt_bar.pdf',
       path = '../plots/controls',
       plot = p_controls_10diff3_1_12_Lparams_pct_total_bar &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 12)


# --- Linearity Of Controls -----------------------------------------------

controls_linearity = read_tsv('controls_linearity_cols1-4.tsv',
                              col_names = c('library', 'barcode', 'copy_num', 'reads'),
                              col_types = 'ccdd',
                              skip = 1) %>%
    mutate(library = str_replace(library, ':', '_')) %>%
    mutate(library = rename_10P[library]) %>%
    left_join(all_ten_plasmids %>% filter(cluster == 'L3')) %>%
    ## left_join(all_ten_plasmids %>% filter(cluster == 'Raw')) %>%
    replace_na(list(count = 0)) %>%
    select(- cluster)

controls_w_gte3_x_values = c('10P-A', '10P-B', '10P-F', '10P-G')

main_text_controls_w_gte3_x_values = c('10P-A', '10P-G')

r2_for_controls_linearity = controls_linearity %>%
    ## filter(library %in% setdiff(controls_w_gte3_x_values, main_text_controls_w_gte3_x_values)
    filter(library %in% main_text_controls_w_gte3_x_values
           & copy_num > 0 & count > 0) %>%
    nest(data = -library) %>%
    mutate(fit = map(data, ~ lm(log10(count) ~ log10(copy_num), data = .x)),
           glanced = map(fit, glance)
           ) %>%
    unnest(glanced) %>%
    select(library, r.squared) %>%
    mutate(r.squared = signif(r.squared, digits = 2)) %>%
    deframe()

controls_linearity %>%
    ## filter(library %in% setdiff(controls_w_gte3_x_values, main_text_controls_w_gte3_x_values)
    filter(library %in% main_text_controls_w_gte3_x_values
           & copy_num > 0 & count > 0) %>%
    ## filter(library %in% c('10P-A', '10P-G') & copy_num > 0 & count > 0) %>%
    ## mutate(library_w_r2 = expression(atop(library, R^2 = r2_for_controls_linearity[library]))) %>%
    mutate(library_w_r2 = paste0(library, "\nR^2 = ", r2_for_controls_linearity[library])) %>%
    ggplot(aes(x = copy_num, y = reads)) +
    geom_smooth(method = 'lm', se = 0, color = 'gray85') +
    geom_point(color = 'cornflowerblue') +
    ## stat_regline_equation(label.x = 1.2, label.y.npc = c(0.95, 0.95, 0.95, 0.95), aes(label = ..rr.label..)) +
    ## facet_wrap(vars(library), scales = 'free', labeller = label_parsed) +
    facet_wrap(vars(library_w_r2), scales = 'free') +
    scale_x_continuous(trans = 'log10',
                       name = 'Plasmid Copies',
                       minor_breaks = NULL,
                       labels = scales::trans_format('log10', scales::math_format(10^.x))) +
    scale_y_continuous(trans = 'log10',
                       name = 'Read Counts',
                       labels = scales::trans_format('log10', scales::math_format(10^.x))) +
    theme_minimal() +
    theme(panel.spacing = unit(1.5, 'lines'),
          strip.background = element_rect(fill = 'gray90', color = 'gray90'),
          ) ->
    ## p_10P_B_F_linearity
    p_10P_A_G_linearity
    ## p_10P_A_B_F_G_linearity

ggsave('linearity_plots_L3_10P-A_10P-G_scatter.pdf',
       path = '../plots/controls',
       plot = p_10P_A_G_linearity,
       height = 3.5)

ggsave('linearity_plots_L3_10P-B_10P-F_scatter.pdf',
       path = '../plots/controls',
       plot = p_10P_B_F_linearity,
       height = 3.5)

ggsave('linearity_plots_10P_scatter.pdf',
       path = '../plots/controls',
       plot = p_10P_A_B_F_G_linearity)


# --- Compare Message-Passing vs Spherical With L = 2 ---------------------

main_text_10P_controls = c('10P-A', '10P-D', '10P-G')

mp_sphere_comparison = expand.grid(library = c('10diff1', '10diff2', '10diff3', '1_01', '1_12', '10_10', '10_12'),
                                   cluster = c('mp', 'sphere'),
                                   L = c('Raw', 'L1', 'L2', 'L3'),
                                   stringsAsFactors = FALSE) %>%
    mutate(data = pmap(list(library, cluster, L),
                       ~ if (..3 != 'Raw') {read_tsv(paste(..1, ..2, ..3, 'clusters.tsv', sep = '_'),
                                  col_names = c('barcode', 'count', 'elements'),
                                  progress = FALSE) %>%
                           select(-elements)
                         }
                         else {read_tsv(paste(..1, 'Raw_clusters.tsv', sep = '_'),
                                        col_names = c('barcode', 'count'),
                                        progress = FALSE)
                             }
                       )
           ) %>%
    unnest(., data) %>%
    mutate(library = rename_10P[library]) %>%
    group_by(library, cluster, L) %>%
    mutate(frac = count / sum(count)) %>%
    ungroup()

## mp_sphere_selected_barcodes = mp_sphere_comparison %>%
##     group_by(library, cluster) %>%
##     mutate(frac = count / sum(count)) %>%
##     ungroup() %>%
##     split(.$library) %>%
##     imap(~ select(.x, - library)) %>%
##     imap(~ filter(.x, ! barcode %in% ten_control_barcodes)) %>%
##     imap(~ group_by(.x, cluster)) %>%
##     imap(~ slice_max(.x, order_by = frac, n = 2)) %>%
##     imap(~ pull(.x, barcode)) %>%
##     imap(~ union(.x, ten_control_barcodes)) %>%
##     enframe(name = 'library', value = 'selected_barcodes') %>%
##     unnest(cols = selected_barcodes)

## mp_sphere_error_selected_barcodes = mp_sphere_comparison %>%
##     group_by(library, cluster, L) %>%
##     filter(! barcode %in% ten_control_barcodes) %>%
##     slice_max(order_by = frac, n = 3) %>%
##     group_by(library, L) %>%
##     distinct(library, L, barcode)

## mp_sphere_selected_barcodes = expand_grid(
##     library = rename_10P,
##     L = c('L1', 'L2', 'L3', 'Raw'),
##     barcode = ten_control_barcodes) %>%
##     bind_rows(mp_sphere_error_selected_barcodes)

mp_sphere_error_selected_barcodes = mp_sphere_comparison %>%
    group_by(library, cluster, L) %>%
    filter(! barcode %in% ten_control_barcodes) %>%
    slice_max(order_by = frac, n = 3) %>%
    ungroup(cluster) %>%
    select(library, L, barcode) %>%
    distinct() %>%
    ungroup() %>%
    expand_grid(cluster = c('mp', 'sphere'))

mp_sphere_control_barcodes = expand_grid(library = rename_10P,
                                         cluster = c('mp', 'sphere'),
                                         L = c('Raw', 'L1', 'L2', 'L3'),
                                         barcode = ten_control_barcodes)

mp_sphere_selected_barcodes = bind_rows(mp_sphere_control_barcodes,
                                        mp_sphere_error_selected_barcodes)

mp_sphere_comparison %>%
    complete(library, cluster, barcode, L, fill = list(count = 0)) %>%
    group_by(library, cluster) %>%
    mutate(frac = count / sum(count)) %>%
    ungroup() %>%
    right_join(mp_sphere_selected_barcodes %>% filter(library %in% main_text_10P_controls),
               by = c('library' = 'library', 'barcode' = 'selected_barcodes')) %>% 
    group_by(library) %>%
    nest() %>%
    mutate(plot = map(data, ~ ggplot(data = .x, aes(x = reorder(barcode, frac), y = frac)) +
                                geom_col(aes(fill = barcode)) +
                                coord_flip() +
                                facet_grid(cols = vars(factor(cluster,
                                                              levels = c('mp', 'sphere'),
                                                              labels = c('Message-Passing', 'Spherical'))),
                                           rows = vars(factor(L, levels = c('L1', 'L2')))) +
                                scale_fill_manual(name = 'Controls',
                                                  values = c(colors_ten_control_barcodes, 'gray85'),
                                                  limits = c(ten_control_barcodes, 'Other'),
                                                  na.value = 'gray85') +
                                scale_y_continuous('Percent Of Total Barcodes',
                                                   labels = scales::label_percent(accuracy = 1),
                                                   ## breaks = c(.125, 0.25, 0.375, 0.5),
                                                   minor_breaks = NULL) +
                                labs(tag = library) +
                                theme_minimal() +
                                theme(panel.spacing.y = unit(2, "lines"),
                                      axis.title = element_blank(),
                                      axis.text.y = element_text(family = 'IBM_Plex_Mono'),
                                      legend.text = element_text(family = 'IBM_Plex_Mono'),
                                      legend.position = 'none'
                                )
                      )
           ) %>%
    pull(plot) %>%
    reduce(`/`) +
    ## plot_annotation(
    ##     tag_levels = 'A',
    ##     tag_prefix = 'Fig. ',
    ##     theme = theme(plot.title = element_text(size = 18))
    ## ) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_ADG_pct_bar
    ## p_10P_A_D_G_L_params_pct_bar
    ## p_10P_C_D_L_params_pct_bar
    ## p_10P_L_params_pct_bar

mp_sphere_plots = mp_sphere_comparison %>%
    semi_join(mp_sphere_selected_barcodes) %>%
    group_by(library) %>%
    nest() %>%
    ## ungroup() %>%
    mutate(plot = map(data, ~ ggplot(data = .x, aes(x = reorder(barcode, frac), y = frac)) +
                                geom_col(aes(fill = barcode)) +
                                coord_flip() +
                                facet_grid(cols = vars(factor(cluster,
                                                              levels = c('mp', 'sphere'),
                                                              labels = c('Message-Passing', 'Spherical'))),
                                           rows = vars(factor(L, levels = c('Raw', 'L1', 'L2', 'L3'))),
                                           scales = 'free') +
                                scale_fill_manual(name = 'Controls',
                                                  values = c(colors_ten_control_barcodes, 'gray85'),
                                                  limits = c(ten_control_barcodes, 'Other'),
                                                  na.value = 'gray85') +
                                scale_y_continuous('Percent of Total Barcodes',
                                                   labels = scales::label_percent(accuracy = 1),
                                                   ## breaks = c(.125, 0.25, 0.375, 0.5),
                                                   minor_breaks = NULL) +
                                ## scale_x_discrete('Barcodes') +
                                labs(title = library,
                                     x = 'Barcodes') +
                                theme_minimal() +
                                theme(panel.spacing.y = unit(2, "lines"),
                                      strip.text = element_text(size = 12),
                                      axis.title = element_blank(),
                                      axis.text.y = element_text(family = 'IBM_Plex_Mono'),
                                      axis.text.x = element_text(size = 12),
                                      legend.text = element_text(family = 'IBM_Plex_Mono'),
                                      legend.position = 'none'
                                )
                      )
           ) %>%
    ungroup() %>%
    select(- data) %>%
    deframe()

## Plot mp vs spherical for 10P-D
mp_sphere_plots %>%
    filter(library == '10P-D') %>%
    pull(plot) %>%
    pluck(1) -> p_mp_sphere_10P_D

p_mp_sphere_10P_D = mp_sphere_plots[['10P-D']]

ggsave('controls_mp_sphere_10P_D_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_D &
           scale_y_continuous('Percent of Total Test',
                              trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              minor_breaks = NULL),
       width = 10,
       height = 17)

ggsave('controls_mp_sphere_10P_ADG_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_ADG_pct_bar &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 10,
       height = 15)

mp_sphere_plots %>%
    ## pull(plot) %>%
    reduce(`/`) +
    ## plot_annotation(
    ##     tag_levels = 'A',
    ##     tag_prefix = 'Fig. ',
    ##     theme = theme(plot.title = element_text(size = 18))
    ## ) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_all_pct_bar

p_10P_mp_sphere_bar = mp_sphere_plots
p_10P_ADG_bar[['10P-A']] = p_10P_ADG_bar[['10P-A']] +
    theme(axis.title.x = element_blank()) #+
    ## theme(axis.title.y = element_blank())
p_10P_ADG_bar[['10P-D']] = p_10P_ADG_bar[['10P-D']] +
    theme(axis.title.x = element_blank())



ggsave('controls_mp_sphere_10P_all_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_all_pct_bar &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL) &
           scale_x_discrete('Barcodes'),
       width = 10,
       height = 40)


## Try ABC plots side-by-side using patchwork
mp_sphere_plots %>%
    slice(4:7) %>%
    pull(plot) %>%
    ## plot_annotation(
    ##     tag_levels = 'A',
    ##     tag_prefix = 'Fig. ',
    ##     theme = theme(plot.title = element_text(size = 18))
    ## ) +
    reduce(`+`) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_DEFG

ggsave('controls_mp_sphere_10P_DEFG_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_DEFG &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 21,
       height = 15)

## Plots for manuscript: ADG for main text, BCEF for supplemental

# ADG
mp_sphere_plots %>%
    filter(library %in% c('10P-A', '10P-D', '10P-G')) %>%
    pull(plot) %>%
    ## plot_annotation(
    ##     tag_levels = 'A',
    ##     tag_prefix = 'Fig. ',
    ##     theme = theme(plot.title = element_text(size = 18))
    ## ) +

mp_sphere_plots[c('10P-A', '10P-D', '10P-G')] %>%
    reduce(`+`) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_ADG

ggsave('controls_mp_sphere_10P_ADG_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_ADG &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 21,
       height = 14)

# BCEF
mp_sphere_plots %>%
    filter(library %in% c('10P-B', '10P-C', '10P-E', '10P-F')) %>%
    pull(plot) %>%
    ## plot_annotation(
    ##     tag_levels = 'A',
    ##     tag_prefix = 'Fig. ',
    ##     theme = theme(plot.title = element_text(size = 18))
    ## ) +
    reduce(`+`) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_BCEF

ggsave('controls_mp_sphere_10P_BCEF_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_BCEF &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 14,
       height = 28)

# DG
mp_sphere_plots %>%
    filter(library %in% c('10P-D', '10P-G')) %>%
    pull(plot) %>%
    reduce(`+`) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_DG

ggsave('controls_mp_sphere_10P_DG_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_DG &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 14)

# ABCEF
mp_sphere_plots %>%
    filter(library %in% c('10P-A', '10P-B', '10P-C', '10P-E', '10P-F')) %>%
    pull(plot) %>%
    reduce(`+`) +
    plot_layout(nrow = 2) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_ABCEF

ggsave('controls_mp_sphere_10P_ABCEF_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_ABCEF &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 22,
       height = 20)


# All Of Them
mp_sphere_plots %>%
    ## pull(plot) %>%
    reduce(`+`) +
    plot_layout(nrow = 2) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_all

ggsave('controls_mp_sphere_10P_all_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_all &
           scale_y_continuous(name = 'Percent of Total Barcodes',
                              trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL) &
           scale_x_discrete('Barcodes'),
       width = 21,
       height = 28)


## Calculate what percentage of barcodes aren't the canonical barcode
one_plasmid_counts %>%
    group_by(library, cluster) %>%
    summarize(non_canonical_pct =
                  sum((barcode != ten_control_barcodes[['BC7_1']]) * count) /
                  sum(count)
              ) %>%
    mutate(non_canonical_pct = scales::percent(non_canonical_pct, accuracy = 0.01)) %>%
    pivot_wider(names_from = cluster, values_from = non_canonical_pct) %>%
    left_join(tibble(library = c('BC1_102', 'BC1_103', 'BC1_104', 'BC1_105'),
                     library_name = c('1P-A', '1P-B', '1P-C', '1P-D'))  # This is to rename libraries
              ) %>%
    ungroup() %>%
    ## arrange(match(cluster, c('Raw', 'L1', 'L2', 'L3')), .by_group = TRUE) %>% # Replaced by pivot_wider, left_joinselect
    select(library_name, Raw, L1, L2, L3) %>%
    kbl(caption = 'Percent Non-Plasmid Barcodes By Levenshtein Parameter',
        col.names = c('Library', 'Raw', 'L1', 'L2', 'L3'),
        align = 'lrrrr') %>%
    ## column_spec(1, bold = TRUE) %>%
    ## collapse_rows(columns = 1, valign = 'top') %>%
    kable_styling(bootstrap_options='striped', full_width=FALSE) %>%
    save_kable('../plots/controls/tables/one_plasmid_pct_non_canonical_bcs.html')


one_plasmid_counts %>%
    filter(cluster == 'Raw') %>%
    group_by(library) %>%
    mutate(frac = count / sum(count)) %>%
    arrange(library, desc(frac)) %>%
    ungroup() %>%
    group_by(barcode) %>%
    mutate(max_frac_across_libraries = max(frac)) %>%
    ungroup() %>%
    filter(barcode %in% (top_n(., 60, wt = max_frac_across_libraries) %>%
                         pull(barcode) %>%
                         unique())
           ) %>%
    arrange(desc(max_frac_across_libraries)) %>%
    mutate(barcode = fct_reorder(as.factor(barcode), max_frac_across_libraries)) %>%
    left_join(tibble(library = c('BC1_102', 'BC1_103', 'BC1_104', 'BC1_105'),
                     library_name = c('1P-A', '1P-B', '1P-C', '1P-D'))  # This is to rename libraries
        ) %>%
    ggplot(aes(x = barcode, y = frac)) +
    geom_col(aes(fill = barcode)) +
    coord_flip() +
    facet_grid(cols = vars(library_name)) +
    ## labs(title = 'One-Plasmid Controls Unclustered', 
    ##      subtitle = 'Fraction Of Total Reads',
    ##      x = 'Barcode') +
    ## scale_fill_manual(name = 'Control',
    ##                   values = colors_ten_control_barcodes,
    ##                   limits = c(ten_control_barcodes[['BC7_1']], 'Other'),
    ##                   na.value = 'gray85') +
    labs(x = 'Barcodes') +
    scale_fill_manual(name = 'Barcode',
                      values = colors_ten_control_barcodes,
                      limits = c(ten_control_barcodes[['BC7_1']], 'Error'),
                      na.value = 'gray85') +
    scale_y_continuous(labels = scales::label_percent(),
                       minor_breaks = NULL) +
    theme_minimal() +
    theme(panel.spacing.y = unit(2, "lines"),
          ## axis.title.x = element_blank(),
          axis.text.y = element_text(family = 'IBM_Plex_Mono'),
          axis.title.y = element_text(vjust = 2),
          axis.title.x = element_text(vjust = -1),
        legend.text = element_text(family = 'IBM_Plex_Mono'),
          ## legend.position = 'none',
          plot.caption = element_text(hjust = 0, size = 12)
          ) -> p_one_plasmid_pct_total_bar

ggsave('one_plasmid_pct_unclustered_bar.pdf',
       path = '../plots/controls',
       plot = p_one_plasmid_pct_total_bar,
       width = 14)
       
ggsave('one_plasmid_pct_unclustered_sqrt_bar.pdf',
       plot = p_one_plasmid_pct_total_bar +
           ## scale_y_continuous('Percent Of Total Barcodes',
           ##                    trans = 'sqrt',
           ##                    labels = scales::label_percent(),
           ##                    minor_breaks = NULL) +
           scale_y_continuous('Square Root-Transformed Percentage of Total Read Counts',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              minor_breaks = NULL),
       path='../plots/controls',
       width=14)

# EPS Version
ggsave('one_plasmid_pct_unclustered_sqrt_bar.eps',
       plot = p_one_plasmid_pct_total_bar +
           ## scale_y_continuous('Percent Of Total Barcodes',
           ##                    trans = 'sqrt',
           ##                    labels = scales::label_percent(),
           ##                    minor_breaks = NULL) +
           scale_y_continuous('Square Root-Transformed Percentage Of Total Read Counts',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              minor_breaks = NULL),
       path='../plots/controls',
       device = 'eps',
       width=14)

# TIFF Version
ggsave('one_plasmid_pct_unclustered_sqrt_bar.tiff',
       plot = p_one_plasmid_pct_total_bar +
           ## scale_y_continuous('Percent Of Total Barcodes',
           ##                    trans = 'sqrt',
           ##                    labels = scales::label_percent(),
           ##                    minor_breaks = NULL) +
           scale_y_continuous('Square Root-Transformed Percentage Of Total Read Counts',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              minor_breaks = NULL),
       path='../plots/controls',
       device = 'tiff',
       dpi = 300,
       compression = 'lzw',
       width = 2200,
       units = 'px')

## -- BC1_103 -------------------------------------------------------------

one_plasmid_counts %>%
    filter(library == 'BC1_103') %>%
    group_by(cluster) %>%
    mutate(frac = count / sum(count)) %>%
    arrange(cluster, desc(frac)) %>%
    ungroup() %>%
    group_by(barcode) %>%
    mutate(max_frac_across_clusters = max(frac)) %>%
    ungroup() %>%
    filter(barcode %in% (group_by(., barcode) %>%
                         filter(frac == max(frac)) %>%
                         ungroup() %>%
                         top_n(15, wt = frac) %>%
                         pull(barcode))
           ) %>%
    complete(library, cluster, nesting(barcode), fill = list(count = 0, frac = 0, max_frac_across_clusters = 0)) %>%
    ## arrange(desc(max_frac_across_clusters)) %>%
    ## mutate(barcode = fct_inorder(as.factor(barcode), ordered = FALSE)) %>%
    ggplot(aes(x = reorder(barcode, max_frac_across_clusters), y = frac)) +
    geom_col(aes(fill = barcode)) +
    coord_flip() +
    facet_grid(cols = vars(factor(cluster, levels = c('Raw', 'L1', 'L2', 'L3')))) +
    labs(title = 'One-Plasmid Control BC1_103', 
         subtitle = 'Fraction Of Total Reads',
         x = 'Barcode') +
    scale_fill_manual(name = 'Control',
                      values = colors_ten_control_barcodes,
                      limits = c(ten_control_barcodes[['BC7_1']], 'Other'),
                      na.value = 'gray85') +
    scale_y_continuous('Percent Of Total Barcodes',
                       labels = scales::label_percent(),
                       minor_breaks = NULL) +
    theme_minimal() +
    theme(panel.spacing.y = unit(2, "lines"),
          axis.text.y = element_text(family = 'IBM_Plex_Mono'),
          axis.title = element_blank(),
          legend.text = element_text(family = 'IBM_Plex_Mono'),
          text = element_text(family = 'Arial')
          ) -> p_one_plasmid_BC1_103_Lparams_pct_total_bar

ggsave('one_plasmid_BC1_103_pct_clustered_bar.pdf',
       path = '../plots/controls',
       plot = p_one_plasmid_BC1_103_Lparams_pct_total_bar,
       width = 14)
       
ggsave('one_plasmid_BC1_103_pct_clustered_sqrt_bar.pdf',
       plot = last_plot() +
           scale_y_continuous('Percent Of Total Barcodes',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              minor_breaks = NULL) +
           labs(subtitle='Fraction Of Total Reads, Square Root'),
       path='../plots/controls',
       width=14)

# EPS/Arial Version
ggsave('one_plasmid_BC1_103_pct_clustered_sqrt_bar.eps',
       plot = last_plot() +
           scale_y_continuous('Percent Of Total Barcodes',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              minor_breaks = NULL) +
           labs(subtitle='Fraction Of Total Reads, Square Root'),
       path='../plots/controls',
       device = 'eps',
       width=14)

# --- 10diff1 Controls ----------------------------------------------------

diff_plasmid_counts = expand.grid(library = c('10diff1', '10diff2', '10diff3'),
                                  cluster = c('Raw', 'L1', 'L2', 'L3'),
                                  stringsAsFactors = FALSE) %>%
    ## filter(cluster != 'Raw') %>%
    mutate(data = map2(.x = library, .y = cluster,
                       ~ if (.y != 'Raw') {
                             read_tsv(paste(.x, 'mp', .y, 'clusters.tsv', sep = '_'),
                                      col_names = c('barcode', 'count', 'elements'),
                                      progress = FALSE) %>%
                                 select(-elements)
                         }
                         else {
                             read.table(paste0('../barcodes_e0.1/JA19375/', .x, '_cutadapt_counts.txt'),
                                              header = FALSE,
                                              col.names = c('count', 'barcode'),
                                              stringsAsFactors = FALSE)
                               })) %>%
    unnest(., data)

### 10diff1 (10P-A)

diff_plasmid_counts %>%
    filter(library == '10diff3') %>%
    group_by(cluster) %>%
    mutate(frac = count / sum(count)) %>%
    arrange(cluster, desc(frac)) %>%
    ungroup() %>%
    group_by(barcode) %>%
    mutate(max_frac_across_clusters = max(frac)) %>%
    ungroup() %>%
    ## filter(barcode %in% (group_by(., barcode) %>%
    ##                      filter(frac == max(frac)) %>%
    ##                      ungroup() %>%
    ##                      top_n(15, wt = frac) %>%
    ##                      pull(barcode))
    ##        ) %>%
    filter(barcode %in% union(as.vector(ten_control_barcodes),
                                        group_by(., barcode) %>%
                                        filter(frac == max(frac)) %>%
                                        ungroup() %>%
                                        filter(!barcode %in% ten_control_barcodes) %>%
                                        top_n(5, wt = frac) %>%
                                        pull(barcode)
                              )
           ) %>%
    complete(library, cluster, nesting(barcode), fill = list(count = 0, frac = 0, max_frac_across_clusters = 0)) %>%
    ggplot(aes(x = reorder(barcode, max_frac_across_clusters), y = frac)) +
    geom_col(aes(fill = barcode)) +
    coord_flip() +
    facet_grid(cols = vars(factor(cluster, levels = c('Raw', 'L1', 'L2', 'L3')))) +
    labs(x = 'Barcode') +
    scale_fill_manual(name = 'Controls',
                      values = c(colors_ten_control_barcodes, 'gray85'),
                      limits = c(ten_control_barcodes, 'Other'),
                      na.value = 'gray85') +
    scale_y_continuous('Percent Of Total Barcodes',
                       labels = scales::label_percent(accuracy = 0.1),
                       breaks = c(0.02, 0.1, 0.25, 0.5),
                       minor_breaks = NULL) +
    theme_minimal() +
    theme(panel.spacing.y = unit(2, "lines"),
          axis.title = element_blank(),
          axis.text.y = element_text(family = 'IBM_Plex_Mono'),
          legend.text = element_text(family = 'IBM_Plex_Mono'),
          legend.position = 'none'
          ) -> p_diff_plasmid_10diff3_Lparams_pct_total_bar

ggsave('diff_plasmid_10diff3_pct_clustered_sqrt_bar.pdf',
       plot = p_diff_plasmid_10diff3_Lparams_pct_total_bar +
           scale_y_continuous('Percent Of Total Barcodes',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              breaks = c(.125, 0.25, 0.375, 0.5),
                              minor_breaks = NULL),
       path='../plots/controls',
       width=14)


diff_plasmid_counts %>%
    filter(library == '10diff2') %>%
    group_by(cluster) %>%
    mutate(frac = count / sum(count)) %>%
    arrange(cluster, desc(frac)) %>%
    ungroup() %>%
    group_by(barcode) %>%
    mutate(max_frac_across_clusters = max(frac)) %>%
    ungroup() %>%
    ## filter(barcode %in% (group_by(., barcode) %>%
    ##                      filter(frac == max(frac)) %>%
    ##                      ungroup() %>%
    ##                      top_n(15, wt = frac) %>%
    ##                      pull(barcode))
    ##        ) %>%
    filter(barcode %in% union(as.vector(ten_control_barcodes),
                                        group_by(., barcode) %>%
                                        filter(frac == max(frac)) %>%
                                        ungroup() %>%
                                        filter(!barcode %in% ten_control_barcodes) %>%
                                        top_n(5, wt = frac) %>%
                                        pull(barcode)
                              )
           ) %>%
    complete(library, cluster, nesting(barcode), fill = list(count = 0, frac = 0, max_frac_across_clusters = 0)) %>%
    ggplot(aes(x = reorder(barcode, max_frac_across_clusters), y = frac)) +
    geom_col(aes(fill = barcode)) +
    coord_flip() +
    facet_grid(cols = vars(factor(cluster, levels = c('Raw', 'L1', 'L2', 'L3')))) +
    labs(x = 'Barcode') +
    scale_fill_manual(name = 'Controls',
                      values = c(colors_ten_control_barcodes, 'gray85'),
                      limits = c(ten_control_barcodes, 'Other'),
                      na.value = 'gray85') +
    scale_y_continuous('Percent Of Total Barcodes',
                       labels = scales::label_percent(accuracy = 0.1),
                       breaks = c(0.02, 0.1, 0.25, 0.5),
                       minor_breaks = NULL) +
    theme_minimal() +
    theme(panel.spacing.y = unit(2, "lines"),
          axis.title = element_blank(),
          axis.text.y = element_text(family = 'IBM_Plex_Mono'),
          legend.text = element_text(family = 'IBM_Plex_Mono'),
          legend.position = 'none'
          ) -> p_diff_plasmid_10diff2_Lparams_pct_total_bar

ggsave('diff_plasmid_10diff2_pct_clustered_bar.pdf',
       path = '../plots/controls',
       plot = p_diff_plasmid_10diff2_Lparams_pct_total_bar,
       width = 14)
       
ggsave('diff_plasmid_10diff2_pct_clustered_sqrt_bar.pdf',
       plot = p_diff_plasmid_10diff2_Lparams_pct_total_bar +
           scale_y_continuous('Percent Of Total Barcodes',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              breaks = c(.125, 0.25, 0.375, 0.5),
                              minor_breaks = NULL),
       path='../plots/controls',
       width=14)


# --- Ten Plasmid Controls ------------------------------------------------

ten_plasmid_counts = expand.grid(library = c('1_01', '1_12', '10_10', '10_12'),
                                 cluster = c('Raw', 'L1', 'L2', 'L3'),
                                 stringsAsFactors = FALSE) %>%
    ## filter(cluster != 'Raw') %>%
    mutate(data = map2(.x = library, .y = cluster,
                       ~ if (.y != 'Raw') {
                             read_tsv(paste(.x, 'mp', .y, 'clusters.tsv', sep = '_'),
                                      col_names = c('barcode', 'count', 'elements'),
                                      progress = FALSE) %>%
                                 select(-elements)
                         }
                         else {
                             read.table(paste0('../barcodes_e0.1/JA19375/', .x, '_cutadapt_counts.txt'),
                                        header = FALSE,
                                        col.names = c('count', 'barcode'),
                                        fill = TRUE,  # for empty barcodes
                                        stringsAsFactors = FALSE)
                               })) %>%
    unnest(., data)

ten_plasmid_counts %>%
    filter(library == '1_12') %>%
    group_by(cluster) %>%
    mutate(frac = count / sum(count)) %>%
    arrange(cluster, desc(frac)) %>%
    ungroup() %>%
    group_by(barcode) %>%
    mutate(max_frac_across_clusters = max(frac)) %>%
    ungroup() %>%
    ## filter(barcode %in% (group_by(., barcode) %>%
    ##                      filter(frac == max(frac)) %>%
    ##                      ungroup() %>%
    ##                      top_n(15, wt = frac) %>%
    ##                      pull(barcode))
    ##        ) %>%
    filter(barcode %in% union(as.vector(ten_control_barcodes),
                                        group_by(., barcode) %>%
                                        filter(frac == max(frac)) %>%
                                        ungroup() %>%
                                        filter(!barcode %in% ten_control_barcodes) %>%
                                        top_n(5, wt = frac) %>%
                                        pull(barcode)
                              )
           ) %>%
    complete(library, cluster, nesting(barcode), fill = list(count = 0, frac = 0, max_frac_across_clusters = 0)) %>%
    ggplot(aes(x = reorder(barcode, max_frac_across_clusters), y = frac)) +
    geom_col(aes(fill = barcode)) +
    coord_flip() +
    facet_grid(cols = vars(factor(cluster, levels = c('Raw', 'L1', 'L2', 'L3')))) +
    ## labs(title = 'Ten??? Control 1_12',
    ##      subtitle = 'Fraction Of Total Reads',
    ##      x = 'Barcode') +
    scale_fill_manual(name = 'Controls',
                      values = c(colors_ten_control_barcodes, 'gray85'),
                      limits = c(ten_control_barcodes, 'Other'),
                      na.value = 'gray85') +
    scale_y_continuous('Percent Of Total Barcodes',
                       labels = scales::label_percent(accuracy = 1),
                       ## breaks = c(.125, 0.25, 0.375, 0.5),
                       minor_breaks = NULL) +
    theme_minimal() +
    theme(panel.spacing.y = unit(2, "lines"),
          axis.title = element_blank(),
          axis.text.y = element_text(family = 'IBM_Plex_Mono'),
          legend.text = element_text(family = 'IBM_Plex_Mono'),
          legend.position = 'none'
          ) -> p_ten_plasmid_1_12_Lparams_pct_total_bar

ggsave('ten_plasmid_1_12_pct_clustered_bar.pdf',
       path = '../plots/controls',
       plot = p_ten_plasmid_1_12_Lparams_pct_total_bar,
       width = 14)
       
ggsave('ten_plasmid_1_12_pct_clustered_sqrt_bar.pdf',
       plot = last_plot() +
           scale_y_continuous('Percent Of Total Barcodes',
                              trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              ## breaks = c(.125, 0.25, 0.375, 0.5),
                              minor_breaks = NULL),
       path='../plots/controls',
       width=14)


# --- Ten-Plasmid Controls All Together -----------------------------------

rename_10P = c(
    '10diff1' = '10P-A',
    '10diff2' = '10P-B',
    '10diff3' = '10P-C',
    '1_01' = '10P-D',
    '1_12' = '10P-E',
    '10_10' = '10P-F',
    '10_12' = '10P-G'
    )

all_ten_plasmids =
    bind_rows(diff_plasmid_counts, ten_plasmid_counts) %>%
    mutate(library = rename_10P[library])

plot_all_ten = all_ten_plasmids %>%
    group_by(library, cluster) %>%
    mutate(frac = count / sum(count)) %>%
    ## arrange(library, cluster, desc(frac)) %>%
    ## ungroup(cluster) %>%
    ## group_by(library, barcode) %>%
    ## ## mutate(max_frac_across_clusters = max(frac)) %>%
    ungroup() %>%
    ## group_by(library) %>%
    ## filter(barcode %in% ten_control_barcodes) %>%
    ## filter(barcode %in% (
    ##                                     group_by(., barcode) %>%
    ##                                     filter(frac == max(frac)) %>%
    ##                                     ungroup() %>%
    ##                                     filter(!barcode %in% ten_control_barcodes) %>%
    ##                                     top_n(5, wt = frac) %>%
    ##                                     pull(barcode)
    ##                           )
    ##        ) %>%
    filter(barcode %in% 
           (
               ## filter(., ! barcode %in% ten_control_barcodes) %>%
     group_by(., library, barcode) %>%
    top_n(1, wt = frac) %>%
    group_by(library) %>%
    top_n(5, wt = frac) %>%
    pull(barcode)
    )
    )
    

    
    filter(barcode %in% union(as.vector(ten_control_barcodes),
                                        group_by(., barcode) %>%
                                        filter(frac == max(frac)) %>%
                                        ungroup() %>%
                                        group_by(library) %>%
                                        filter(!barcode %in% ten_control_barcodes) %>%
                                        top_n(5, wt = frac) %>%
                                        pull(barcode) 
                              )
           )

## attempt with nesting from the beginning

## selected_barcodes = all_ten_plasmids %>%
##     group_by(library, cluster) %>%
##     mutate(frac = count / sum(count)) %>%
##     ungroup() %>%
##     split(.$library) %>%
##     imap(~ select(.x, - library)) %>%
##     imap(~ filter(.x, ! barcode %in% ten_control_barcodes)) %>%
##     imap(~ group_by(.x, barcode)) %>%
##     imap(~ mutate(.x, selection_metric = max(frac))) %>%
##     ## imap(~ mutate(.x, selection_metric = median(frac))) %>%
##     ## imap(~ mutate(.x, selection_metric = mean(frac))) %>%
##     imap(~ distinct(.x, barcode, .keep_all = TRUE)) %>%
##     imap(~ ungroup(.x)) %>%
##     imap(~ top_n(.x, 5, wt = selection_metric)) %>%
##     imap(~ pull(.x, barcode)) %>%
##     imap(~ union(.x, ten_control_barcodes)) %>%
##     enframe(name = 'library', value = 'selected_barcodes') %>%
##     unnest(cols = selected_barcodes)

selected_barcodes = all_ten_plasmids %>%
    group_by(library, cluster) %>%
    mutate(frac = count / sum(count)) %>%
    ungroup() %>%
    split(.$library) %>%
    imap(~ select(.x, - library)) %>%
    imap(~ filter(.x, ! barcode %in% ten_control_barcodes)) %>%
    imap(~ group_by(.x, cluster)) %>%
    imap(~ top_n(.x, 2, wt = frac)) %>%
    imap(~ pull(.x, barcode)) %>%
    imap(~ union(.x, ten_control_barcodes)) %>%
    enframe(name = 'library', value = 'selected_barcodes') %>%
    unnest(cols = selected_barcodes)

main_text_10P_controls = c('10P-A', '10P-D', '10P-G')

all_ten_plasmids %>%
    complete(library, cluster, barcode, fill = list(count = 0)) %>%
    group_by(library, cluster) %>%
    mutate(frac = count / sum(count)) %>%
    ungroup() %>%
    right_join(selected_barcodes %>% filter(library %in% main_text_10P_controls),
               by = c('library' = 'library', 'barcode' = 'selected_barcodes')) %>% 
    group_by(library) %>%
    nest() %>%
    mutate(plot = map(data, ~ ggplot(data = .x, aes(x = reorder(barcode, frac), y = frac)) +
                                geom_col(aes(fill = barcode)) +
                                coord_flip() +
                                facet_grid(cols = vars(factor(cluster, levels = c('Raw', 'L1', 'L2', 'L3')))) +
                                scale_fill_manual(name = 'Controls',
                                                  values = c(colors_ten_control_barcodes, 'gray85'),
                                                  limits = c(ten_control_barcodes, 'Other'),
                                                  na.value = 'gray85') +
                                scale_y_continuous('Percent Of Total Barcodes',
                                                   labels = scales::label_percent(accuracy = 1),
                                                   ## breaks = c(.125, 0.25, 0.375, 0.5),
                                                   minor_breaks = NULL) +
                                labs(tag = library) +
                                theme_minimal() +
                                theme(panel.spacing.y = unit(2, "lines"),
                                      axis.title = element_blank(),
                                      axis.text.y = element_text(family = 'IBM_Plex_Mono'),
                                      legend.text = element_text(family = 'IBM_Plex_Mono'),
                                      legend.position = 'none'
                                )
                      )
           ) %>%
    pull(plot) %>%
    reduce(`/`) +
    ## plot_annotation(
    ##     tag_levels = 'A',
    ##     tag_prefix = 'Fig. ',
    ##     theme = theme(plot.title = element_text(size = 18))
    ## ) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    ## xlab(label = 'Percent Of Total Barcodes') ->
    p_10P_ADG_L_params_pct_bar
    ## p_10P_A_D_G_L_params_pct_bar
    ## p_10P_C_D_L_params_pct_bar
    ## p_10P_L_params_pct_bar

ten_plasmid_L_plots = all_ten_plasmids %>%
    complete(library, cluster, barcode, fill = list(count = 0)) %>%
    group_by(library, cluster) %>%
    mutate(frac = count / sum(count)) %>%
    ungroup() %>%
    semi_join(selected_barcodes, by = c('library' = 'library', 'barcode' = 'selected_barcodes')) %>% 
    group_by(library) %>%
    nest() %>%
    mutate(plot = map(data, ~ ggplot(data = .x, aes(x = reorder(barcode, frac), y = frac)) +
                                geom_col(aes(fill = barcode)) +
                                coord_flip() +
                                facet_grid(cols = vars(factor(cluster, levels = c('Raw', 'L1', 'L2', 'L3')))) +
                                scale_fill_manual(name = 'Controls',
                                                  values = c(colors_ten_control_barcodes, 'gray85'),
                                                  limits = c(ten_control_barcodes, 'Other'),
                                                  na.value = 'gray85') +
                                scale_y_continuous('Percent of Total Barcodes',
                                                   labels = scales::label_percent(accuracy = 1),
                                                   ## breaks = c(.125, 0.25, 0.375, 0.5),
                                                   minor_breaks = NULL) +
                                labs(tag = library) +
                                theme_minimal() +
                                theme(panel.spacing.y = unit(2, "lines"),
                                      ## axis.title = element_blank(),
                                      ## axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12),
          ## axis.title.y = element_text(vjust = 2, size = 12),
          ## axis.title.x = element_text(vjust = -1, size = 12),
          strip.text = element_text(size = 12),
                                      axis.title.y = element_text(vjust = 10, size = 12),
                                      axis.text.y = element_text(family = 'IBM_Plex_Mono'),
                                      axis.title.x = element_text(vjust = -1, size = 12),
                                      legend.text = element_text(family = 'IBM_Plex_Mono'),
                                      legend.position = 'none'
                                )
                      )
           ) %>%
    ungroup() %>%
    select(-data) %>%
    deframe()

p_10P_ADG_bar = ten_plasmid_L_plots[c('10P-A', '10P-D', '10P-G')]
p_10P_ADG_bar[['10P-A']] = p_10P_ADG_bar[['10P-A']] +
    theme(axis.title.x = element_blank()) #+
    ## theme(axis.title.y = element_blank())
p_10P_ADG_bar[['10P-D']] = p_10P_ADG_bar[['10P-D']] +
    theme(axis.title.x = element_blank())
## p_10P_ADG_bar[['10P-G']] = p_10P_ADG_bar[['10P-G']] +
##     theme(axis.title.y = element_blank())

p_10P_ADG_bar %<>%
    reduce(`/`) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12))

ggsave('controls_10P_ADG_bar.pdf',
       path = '../plots/controls',
       plot = p_10P_ADG_bar &
           scale_y_continuous('Square Root-Transformed Percentage of Total Read Counts',
                              trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL) &
           scale_x_discrete('Barcodes'),
       width = 15,
       height = 15
       )

p_10P_BCEF_bar = ten_plasmid_L_plots[c('10P-B', '10P-C', '10P-E', '10P-F')]
p_10P_BCEF_bar[['10P-B']] = p_10P_BCEF_bar[['10P-B']] +
    theme(axis.title.x = element_blank()) #+
    ## theme(axis.title.y = element_blank())
p_10P_BCEF_bar[['10P-C']] = p_10P_BCEF_bar[['10P-C']] +
    theme(axis.title.x = element_blank())
p_10P_BCEF_bar[['10P-E']] = p_10P_BCEF_bar[['10P-E']] +
    theme(axis.title.x = element_blank())
## p_10P_BCEF_bar[['10P-E']] = p_10P_BCEF_bar[['10P-E']] +
##     theme(axis.title.y = element_blank())

p_10P_BCEF_bar %<>%
    reduce(`/`) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12))

ggsave('controls_10P_BCEF_bar.pdf',
       path = '../plots/controls',
       plot = p_10P_BCEF_bar &
           scale_y_continuous('Square Root-Transformed Percentage of Total Read Counts',
                              trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL) &
           scale_x_discrete('Barcodes'),
       width = 15,
       height = 15
       )


ggsave('controls_all_10P_bar.pdf',
       path = '../plots/controls',
       plot = p_10P_L_params_pct_bar &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 30)

ggsave('controls_10P-C_10P-D_bar.pdf',
       path = '../plots/controls',
       plot = p_10P_C_D_L_params_pct_bar &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 10)

ggsave('controls_10P_ADG_bar.pdf',
       path = '../plots/controls',
       plot = p_10P_ADG_L_params_pct_bar +
           xlab(label = 'Percent Of Total Barcodes') &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 15)

ggsave('controls_10P-B-C-E-F_bar.pdf',
       path = '../plots/controls',
       plot = p_10P_B_C_E_F_L_params_pct_bar &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 20)



# --- Compose Single Sequence Of L Of Three Types Of Controls -------------

p_controls_BC1_103_10diff2_1_12_Lparams_pct_total_bar =
    (
        p_one_plasmid_BC1_103_Lparams_pct_total_bar +
        theme(legend.position = 'None')
    ) /
        p_diff_plasmid_10diff2_Lparams_pct_total_bar / 
        p_ten_plasmid_1_12_Lparams_pct_total_bar +
    plot_annotation(
        ## title = 'Percent Of Total Reads For Controls',
        ## subtitle = 'Square Root Transformed',
        tag_levels = 'A',
        tag_prefix = 'Fig. ',
        theme = theme(plot.title = element_text(size = 18))
        ) +
    plot_layout(guides = 'collect') &
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12))

p_controls_10diff2_1_12_Lparams_pct_total_bar =
        p_diff_plasmid_10diff2_Lparams_pct_total_bar / 
        p_ten_plasmid_1_12_Lparams_pct_total_bar +
    plot_annotation(
        ## title = 'Percent Of Total Reads For Controls',
        ## subtitle = 'Square Root Transformed',
        tag_levels = 'A',
        tag_prefix = 'Fig. ',
        theme = theme(plot.title = element_text(size = 18))
        ) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12))

p_controls_10diff3_1_12_Lparams_pct_total_bar =
        p_diff_plasmid_10diff3_Lparams_pct_total_bar / 
        p_ten_plasmid_1_12_Lparams_pct_total_bar +
    plot_annotation(
        ## title = 'Percent Of Total Reads For Controls',
        ## subtitle = 'Square Root Transformed',
        tag_levels = 'A',
        tag_prefix = 'Fig. ',
        theme = theme(plot.title = element_text(size = 18))
        ) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12))

ggsave('controls_BC1_103_10diff2_1_12_clustered_sqrt_bar.pdf',
       path = '../plots/controls',
       plot = p_controls_BC1_103_10diff2_1_12_Lparams_pct_total_bar &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 18)


ggsave('controls_10diff3_1_12_clustered_sqrt_bar.pdf',
       path = '../plots/controls',
       plot = p_controls_10diff3_1_12_Lparams_pct_total_bar &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 12)


# --- Linearity Of Controls -----------------------------------------------

controls_linearity = read_tsv('controls_linearity_cols1-4.tsv',
                              col_names = c('library', 'barcode', 'copy_num', 'reads'),
                              col_types = 'ccdd',
                              skip = 1) %>%
    mutate(library = str_replace(library, ':', '_')) %>%
    mutate(library = rename_10P[library]) %>%
    left_join(all_ten_plasmids %>% filter(cluster == 'L3')) %>%
    ## left_join(all_ten_plasmids %>% filter(cluster == 'Raw')) %>%
    replace_na(list(count = 0)) %>%
    select(- cluster)

controls_w_gte3_x_values = c('10P-A', '10P-B', '10P-F', '10P-G')

main_text_controls_w_gte3_x_values = c('10P-A', '10P-G')

r2_for_controls_linearity = controls_linearity %>%
    ## filter(library %in% setdiff(controls_w_gte3_x_values, main_text_controls_w_gte3_x_values)
    filter(library %in% main_text_controls_w_gte3_x_values
           & copy_num > 0 & count > 0) %>%
    nest(data = -library) %>%
    mutate(fit = map(data, ~ lm(log10(count) ~ log10(copy_num), data = .x)),
           glanced = map(fit, glance)
           ) %>%
    unnest(glanced) %>%
    select(library, r.squared) %>%
    mutate(r.squared = signif(r.squared, digits = 2)) %>%
    deframe()

controls_linearity %>%
    ## filter(library %in% setdiff(controls_w_gte3_x_values, main_text_controls_w_gte3_x_values)
    filter(library %in% main_text_controls_w_gte3_x_values
           & copy_num > 0 & count > 0) %>%
    ## filter(library %in% c('10P-A', '10P-G') & copy_num > 0 & count > 0) %>%
    ## mutate(library_w_r2 = expression(atop(library, R^2 = r2_for_controls_linearity[library]))) %>%
    mutate(library_w_r2 = paste0(library, "\nR^2 = ", r2_for_controls_linearity[library])) %>%
    ggplot(aes(x = copy_num, y = reads)) +
    geom_smooth(method = 'lm', se = 0, color = 'gray85') +
    geom_point(color = 'cornflowerblue') +
    ## stat_regline_equation(label.x = 1.2, label.y.npc = c(0.95, 0.95, 0.95, 0.95), aes(label = ..rr.label..)) +
    ## facet_wrap(vars(library), scales = 'free', labeller = label_parsed) +
    facet_wrap(vars(library_w_r2), scales = 'free') +
    scale_x_continuous(trans = 'log10',
                       name = 'Plasmid Copies',
                       minor_breaks = NULL,
                       labels = scales::trans_format('log10', scales::math_format(10^.x))) +
    scale_y_continuous(trans = 'log10',
                       name = 'Read Counts',
                       labels = scales::trans_format('log10', scales::math_format(10^.x))) +
    theme_minimal() +
    theme(panel.spacing = unit(1.5, 'lines'),
          strip.background = element_rect(fill = 'gray90', color = 'gray90'),
          ) ->
    ## p_10P_B_F_linearity
    p_10P_A_G_linearity
    ## p_10P_A_B_F_G_linearity

ggsave('linearity_plots_L3_10P-A_10P-G_scatter.pdf',
       path = '../plots/controls',
       plot = p_10P_A_G_linearity,
       height = 3.5)

ggsave('linearity_plots_L3_10P-B_10P-F_scatter.pdf',
       path = '../plots/controls',
       plot = p_10P_B_F_linearity,
       height = 3.5)

ggsave('linearity_plots_10P_scatter.pdf',
       path = '../plots/controls',
       plot = p_10P_A_B_F_G_linearity)


# --- Compare Message-Passing vs Spherical With L = 2 ---------------------

main_text_10P_controls = c('10P-A', '10P-D', '10P-G')

mp_sphere_comparison = expand.grid(library = c('10diff1', '10diff2', '10diff3', '1_01', '1_12', '10_10', '10_12'),
                                   cluster = c('mp', 'sphere'),
                                   L = c('Raw', 'L1', 'L2', 'L3'),
                                   stringsAsFactors = FALSE) %>%
    mutate(data = pmap(list(library, cluster, L),
                       ~ if (..3 != 'Raw') {read_tsv(paste(..1, ..2, ..3, 'clusters.tsv', sep = '_'),
                                  col_names = c('barcode', 'count', 'elements'),
                                  progress = FALSE) %>%
                           select(-elements)
                         }
                         else {read_tsv(paste(..1, 'Raw_clusters.tsv', sep = '_'),
                                        col_names = c('barcode', 'count'),
                                        progress = FALSE)
                             }
                       )
           ) %>%
    unnest(., data) %>%
    mutate(library = rename_10P[library]) %>%
    group_by(library, cluster, L) %>%
    mutate(frac = count / sum(count)) %>%
    ungroup()

## mp_sphere_selected_barcodes = mp_sphere_comparison %>%
##     group_by(library, cluster) %>%
##     mutate(frac = count / sum(count)) %>%
##     ungroup() %>%
##     split(.$library) %>%
##     imap(~ select(.x, - library)) %>%
##     imap(~ filter(.x, ! barcode %in% ten_control_barcodes)) %>%
##     imap(~ group_by(.x, cluster)) %>%
##     imap(~ slice_max(.x, order_by = frac, n = 2)) %>%
##     imap(~ pull(.x, barcode)) %>%
##     imap(~ union(.x, ten_control_barcodes)) %>%
##     enframe(name = 'library', value = 'selected_barcodes') %>%
##     unnest(cols = selected_barcodes)

## mp_sphere_error_selected_barcodes = mp_sphere_comparison %>%
##     group_by(library, cluster, L) %>%
##     filter(! barcode %in% ten_control_barcodes) %>%
##     slice_max(order_by = frac, n = 3) %>%
##     group_by(library, L) %>%
##     distinct(library, L, barcode)

## mp_sphere_selected_barcodes = expand_grid(
##     library = rename_10P,
##     L = c('L1', 'L2', 'L3', 'Raw'),
##     barcode = ten_control_barcodes) %>%
##     bind_rows(mp_sphere_error_selected_barcodes)

mp_sphere_error_selected_barcodes = mp_sphere_comparison %>%
    group_by(library, cluster, L) %>%
    filter(! barcode %in% ten_control_barcodes) %>%
    slice_max(order_by = frac, n = 3) %>%
    ungroup(cluster) %>%
    select(library, L, barcode) %>%
    distinct() %>%
    ungroup() %>%
    expand_grid(cluster = c('mp', 'sphere'))

mp_sphere_control_barcodes = expand_grid(library = rename_10P,
                                         cluster = c('mp', 'sphere'),
                                         L = c('Raw', 'L1', 'L2', 'L3'),
                                         barcode = ten_control_barcodes)

mp_sphere_selected_barcodes = bind_rows(mp_sphere_control_barcodes,
                                        mp_sphere_error_selected_barcodes)

mp_sphere_comparison %>%
    complete(library, cluster, barcode, L, fill = list(count = 0)) %>%
    group_by(library, cluster) %>%
    mutate(frac = count / sum(count)) %>%
    ungroup() %>%
    right_join(mp_sphere_selected_barcodes %>% filter(library %in% main_text_10P_controls),
               by = c('library' = 'library', 'barcode' = 'selected_barcodes')) %>% 
    group_by(library) %>%
    nest() %>%
    mutate(plot = map(data, ~ ggplot(data = .x, aes(x = reorder(barcode, frac), y = frac)) +
                                geom_col(aes(fill = barcode)) +
                                coord_flip() +
                                facet_grid(cols = vars(factor(cluster,
                                                              levels = c('mp', 'sphere'),
                                                              labels = c('Message-Passing', 'Spherical'))),
                                           rows = vars(factor(L, levels = c('L1', 'L2')))) +
                                scale_fill_manual(name = 'Controls',
                                                  values = c(colors_ten_control_barcodes, 'gray85'),
                                                  limits = c(ten_control_barcodes, 'Other'),
                                                  na.value = 'gray85') +
                                scale_y_continuous('Percent Of Total Barcodes',
                                                   labels = scales::label_percent(accuracy = 1),
                                                   ## breaks = c(.125, 0.25, 0.375, 0.5),
                                                   minor_breaks = NULL) +
                                labs(tag = library) +
                                theme_minimal() +
                                theme(panel.spacing.y = unit(2, "lines"),
                                      axis.title = element_blank(),
                                      axis.text.y = element_text(family = 'IBM_Plex_Mono'),
                                      legend.text = element_text(family = 'IBM_Plex_Mono'),
                                      legend.position = 'none'
                                )
                      )
           ) %>%
    pull(plot) %>%
    reduce(`/`) +
    ## plot_annotation(
    ##     tag_levels = 'A',
    ##     tag_prefix = 'Fig. ',
    ##     theme = theme(plot.title = element_text(size = 18))
    ## ) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_ADG_pct_bar
    ## p_10P_A_D_G_L_params_pct_bar
    ## p_10P_C_D_L_params_pct_bar
    ## p_10P_L_params_pct_bar

mp_sphere_plots = mp_sphere_comparison %>%
    semi_join(mp_sphere_selected_barcodes) %>%
    group_by(library) %>%
    nest() %>%
    ## ungroup() %>%
    mutate(plot = map(data, ~ ggplot(data = .x, aes(x = reorder(barcode, frac), y = frac)) +
                                geom_col(aes(fill = barcode)) +
                                coord_flip() +
                                facet_grid(cols = vars(factor(cluster,
                                                              levels = c('mp', 'sphere'),
                                                              labels = c('Message-Passing', 'Spherical'))),
                                           rows = vars(factor(L, levels = c('Raw', 'L1', 'L2', 'L3'))),
                                           scales = 'free') +
                                scale_fill_manual(name = 'Controls',
                                                  values = c(colors_ten_control_barcodes, 'gray85'),
                                                  limits = c(ten_control_barcodes, 'Other'),
                                                  na.value = 'gray85') +
                                scale_y_continuous('Percent of Total Barcodes',
                                                   labels = scales::label_percent(accuracy = 1),
                                                   ## breaks = c(.125, 0.25, 0.375, 0.5),
                                                   minor_breaks = NULL) +
                                ## scale_x_discrete('Barcodes') +
                                labs(title = library,
                                     x = 'Barcodes') +
                                theme_minimal() +
                                theme(panel.spacing.y = unit(2, "lines"),
                                      axis.title = element_blank(),
                                      axis.text.y = element_text(family = 'IBM_Plex_Mono'),
                                      legend.text = element_text(family = 'IBM_Plex_Mono'),
                                      legend.position = 'none'
                                )
                      )
           ) %>%
    ungroup() %>%
    select(- data) %>%
    deframe()

## Plot mp vs spherical for 10P-D
mp_sphere_plots %>%
    filter(library == '10P-D') %>%
    pull(plot) %>%
    pluck(1) -> p_mp_sphere_10P_D

p_mp_sphere_10P_D = mp_sphere_plots[['10P-D']]

ggsave('controls_mp_sphere_10P_D_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_D &
           scale_y_continuous('Percent of Total Test',
                              trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              minor_breaks = NULL),
       width = 10,
       height = 17)

ggsave('controls_mp_sphere_10P_ADG_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_ADG_pct_bar &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 10,
       height = 15)

mp_sphere_plots %>%
    pull(plot) %>%
    reduce(`/`) +
    ## plot_annotation(
    ##     tag_levels = 'A',
    ##     tag_prefix = 'Fig. ',
    ##     theme = theme(plot.title = element_text(size = 18))
    ## ) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_all_pct_bar

p_10P_mp_sphere_bar = mp_sphere_plots
p_10P_ADG_bar[['10P-A']] = p_10P_ADG_bar[['10P-A']] +
    theme(axis.title.x = element_blank()) #+
    ## theme(axis.title.y = element_blank())
p_10P_ADG_bar[['10P-D']] = p_10P_ADG_bar[['10P-D']] +
    theme(axis.title.x = element_blank())



ggsave('controls_mp_sphere_10P_all_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_all_pct_bar &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL) &
           scale_x_discrete('Barcodes'),
       width = 10,
       height = 40)


## Try ABC plots side-by-side using patchwork
mp_sphere_plots %>%
    slice(4:7) %>%
    pull(plot) %>%
    ## plot_annotation(
    ##     tag_levels = 'A',
    ##     tag_prefix = 'Fig. ',
    ##     theme = theme(plot.title = element_text(size = 18))
    ## ) +
    reduce(`+`) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_DEFG

ggsave('controls_mp_sphere_10P_DEFG_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_DEFG &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 21,
       height = 15)

## Plots for manuscript: ADG for main text, BCEF for supplemental

# ADG
mp_sphere_plots %>%
    filter(library %in% c('10P-A', '10P-D', '10P-G')) %>%
    pull(plot) %>%
    ## plot_annotation(
    ##     tag_levels = 'A',
    ##     tag_prefix = 'Fig. ',
    ##     theme = theme(plot.title = element_text(size = 18))
    ## ) +
    reduce(`+`) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_ADG

ggsave('controls_mp_sphere_10P_ADG_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_ADG &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 21,
       height = 14)

# BCEF
mp_sphere_plots %>%
    filter(library %in% c('10P-B', '10P-C', '10P-E', '10P-F')) %>%
    pull(plot) %>%
    ## plot_annotation(
    ##     tag_levels = 'A',
    ##     tag_prefix = 'Fig. ',
    ##     theme = theme(plot.title = element_text(size = 18))
    ## ) +
    reduce(`+`) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_BCEF

ggsave('controls_mp_sphere_10P_BCEF_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_BCEF &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 14,
       height = 28)

# DG
mp_sphere_plots %>%
    filter(library %in% c('10P-D', '10P-G')) %>%
    pull(plot) %>%
    reduce(`+`) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_DG

ggsave('controls_mp_sphere_10P_DG_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_DG &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 14)

# ABCEF
mp_sphere_plots %>%
    filter(library %in% c('10P-A', '10P-B', '10P-C', '10P-E', '10P-F')) %>%
    pull(plot) %>%
    reduce(`+`) +
    plot_layout(nrow = 2) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_ABCEF

ggsave('controls_mp_sphere_10P_ABCEF_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_ABCEF &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 22,
       height = 20)


# All Of Them
mp_sphere_plots %>%
    ## pull(plot) %>%
    reduce(`+`) +
    plot_layout(nrow = 2) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_mp_sphere_10P_all

ggsave('controls_mp_sphere_10P_all_bar.pdf',
       path = '../plots/controls',
       plot = p_mp_sphere_10P_all &
           scale_y_continuous(name = 'Percent of Total Barcodes',
                              trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL) &
           scale_x_discrete('Barcodes'),
       width = 21,
       height = 28)


# --- Experiments With Throwing Out "Low-Quality Reads" -------------------

## -- BC1 Plasmids Throwing Out Low-Quality Reads -------------------------

gte20_one_plasmid_counts = tibble(library = c('BC1_102', 'BC1_103', 'BC1_104', 'BC1_105')) %>%
    mutate(barcode_counts = map(library,
                                ~ read_table(paste0('../barcodes_e0.1/cutadapt_fastqs/JA19375/', ., '_qual_gte7_counts.txt'),
                                             col_names = c('count', 'barcode')
                                             )
                                )
           ) %>%
    unnest(barcode_counts)

gte20_one_plasmid_counts %>%
    group_by(library) %>%
    mutate(frac = count / sum(count)) %>%
    arrange(library, desc(frac)) %>%
    ungroup() %>%
    group_by(barcode) %>%
    mutate(max_frac_across_libraries = max(frac)) %>%
    ungroup() %>%
    filter(barcode %in% (top_n(., 60, wt = max_frac_across_libraries) %>%
                         pull(barcode) %>%
                         unique())
           ) %>%
    arrange(desc(max_frac_across_libraries)) %>%
    mutate(barcode = fct_reorder(as.factor(barcode), max_frac_across_libraries)) %>%
    left_join(tibble(library = c('BC1_102', 'BC1_103', 'BC1_104', 'BC1_105'),
                     library_name = c('1P-A', '1P-B', '1P-C', '1P-D'))  # This is to rename libraries
        ) %>%
    ggplot(aes(x = barcode, y = frac)) +
    geom_col(aes(fill = barcode)) +
    coord_flip() +
    facet_grid(cols = vars(library_name)) +
    ## labs(title = 'One-Plasmid Controls Unclustered', 
    ##      subtitle = 'Fraction Of Total Reads',
    ##      x = 'Barcode') +
    ## scale_fill_manual(name = 'Control',
    ##                   values = colors_ten_control_barcodes,
    ##                   limits = c(ten_control_barcodes[['BC7_1']], 'Other'),
    ##                   na.value = 'gray85') +
    labs(x = 'Barcodes') +
    scale_fill_manual(name = 'Barcode',
                      values = colors_ten_control_barcodes,
                      limits = c(ten_control_barcodes[['BC7_1']], 'Error'),
                      na.value = 'gray85') +
    scale_y_continuous(labels = scales::label_percent(),
                       minor_breaks = NULL) +
    theme_minimal() +
    theme(panel.spacing.y = unit(2, "lines"),
          ## axis.title.x = element_blank(),
          axis.text.y = element_text(family = 'IBM_Plex_Mono'),
          axis.title.y = element_text(vjust = 2),
          axis.title.x = element_text(vjust = -1),
        legend.text = element_text(family = 'IBM_Plex_Mono'),
          ## legend.position = 'none',
          plot.caption = element_text(hjust = 0, size = 12)
          ) -> p_gte20_one_plasmid_pct_total_bar

ggsave('one_plasmid_pct_unclustered_bar.pdf',
       path = '../plots/controls',
       plot = p_one_plasmid_pct_total_bar,
       width = 14)
       
ggsave('gte20_one_plasmid_pct_unclustered_sqrt_bar.pdf',
       plot = p_gte20_one_plasmid_pct_total_bar +
           ## scale_y_continuous('Percent Of Total Barcodes',
           ##                    trans = 'sqrt',
           ##                    labels = scales::label_percent(),
           ##                    minor_breaks = NULL) +
           scale_y_continuous('Square Root-Transformed Percentage of Total Read Counts',
                              trans = 'sqrt',
                              labels = scales::label_percent(),
                              minor_breaks = NULL),
       path='../plots/controls',
       width=14)

## Calculate what percentage of barcodes aren't the canonical barcode
gte20_one_plasmid_counts %>%
    group_by(library) %>%
    summarize(non_canonical_pct =
                  sum((barcode != ten_control_barcodes[['BC7_1']]) * count) /
                  sum(count)
              ) %>%
    mutate(non_canonical_pct = scales::percent(non_canonical_pct, accuracy = 0.01)) %>%
    ## arrange(match(cluster, c('Raw', 'L1', 'L2', 'L3')), .by_group = TRUE) %>%
    kbl(caption = 'Percentage Of Non-Canonical Barcodes With All Base PHRED Scores >= 20') %>%
    column_spec(1, bold = TRUE) %>%
    collapse_rows(columns = 1, valign = 'top') %>%
    kable_styling(bootstrap_options='striped', full_width=FALSE) %>%
    save_kable('../plots/controls/tables/gte20_one_plasmid_pct_non_canonical_bcs.html')


## -- 1-Plasmid Controls: Noise Between Raw vs PHRED = 20, 30 --------------

vary_qual_one_plasmid_counts = expand_grid(library = c('BC1_102', 'BC1_103', 'BC1_104', 'BC1_105'),
                                           quality = c('Raw', '20', '30')) %>%
    ## filter(cluster != 'raw') %>%
    mutate(data = map2(.x = library, .y = quality,
                       ~ if (.y != 'Raw') {
                             read_table(paste0('../barcodes_e0.1/cutadapt_fastqs/JA19375/',
                                               .x, '_qual_gte', .y, '_counts.txt'),
                                        col_names = c('count', 'barcode'),
                                        progress = FALSE)
                         }
                         else {
                             read_table(paste0('../barcodes_e0.1/JA19375/',
                                               .x, '_cutadapt_counts.txt'),
                                        col_names = c('count', 'barcode'),
                                        progress = FALSE)})) %>%
    unnest(., data)

vary_qual_one_plasmid_counts %>%
    group_by(library, quality) %>%
    summarize(non_canonical_pct =
                  sum((barcode != ten_control_barcodes[['BC7_1']]) * count) /
                  sum(count)) %>%
    ungroup() %>%
    mutate(non_canonical_pct = scales::percent(non_canonical_pct, accuracy = 0.01)) %>% 
    pivot_wider(names_from = quality, values_from = non_canonical_pct) %>%
    left_join(tibble(library = c('BC1_102', 'BC1_103', 'BC1_104', 'BC1_105'),
                     library_name = c('1P-A', '1P-B', '1P-C', '1P-D'))  # This is to rename libraries
              ) %>%
    select('library_name', 'Raw', '20', '30') %>%
    kbl(caption = 'Percent Non-Plasmid Barcodes With PHRED Score Quality Cutoff',
        col.names = c('Library', 'Raw', '20', '30'),
        align = 'lrrr') %>%
    column_spec(1, bold = TRUE) %>%
    ## collapse_rows(columns = 1, valign = 'top') %>%
    kable_styling(bootstrap_options='striped', full_width=FALSE) %>%
    ## add_header_above(c(' ' = 1, 'Raw' = 2, '20' = 2, '30' = 2)) %>%
    save_kable('../plots/controls/tables/vary_qual_one_plasmid_pct_non_canonical_bcs.html')

    
## Commented version includes number of error barcodes in addition to percent
    ## mutate(non_canonical_pct = scales::percent(non_canonical_pct, accuracy = 0.01),
    ##        num_barcodes = scales::comma(num_barcodes)) %>% 
    ## pivot_wider(names_from = quality,
    ##             values_from = c(non_canonical_pct, num_barcodes)) %>%
    ## left_join(tibble(library = c('BC1_102', 'BC1_103', 'BC1_104', 'BC1_105'),
    ##                  library_name = c('1P-A', '1P-B', '1P-C', '1P-D'))  # This is to rename libraries
    ##     ) %>%
    ## select('library_name', 'non_canonical_pct_Raw', 'num_barcodes_Raw',
    ##        'non_canonical_pct_20', 'num_barcodes_20',
    ##        'non_canonical_pct_30', 'num_barcodes_30') %>%
    ## #
# select('Raw', `20`, `30`) %>%
    ## ## left_join(tibble(library = c('BC1_102', 'BC1_103', 'BC1_104', 'BC1_105'),
    ## ##                  library_name = c('1P-A', '1P-B', '1P-C', '1P-D'))  # This is to rename libraries
    ## ##     ) %>%
    ## kbl(caption = 'Percent Non-Plasmid Barcodes and Num Total Barcodes By Quality Cutoff',
    ##     col.names = c('Library', '% error', '# barcodes', '% error', '# barcodes', '% error', '# barcodes'),
    ##     align = 'lrrrrrr') %>%
    ## column_spec(1, bold = TRUE) %>%
    ## ## collapse_rows(columns = 1, valign = 'top') %>%
    ## kable_styling(bootstrap_options='striped', full_width=FALSE) %>%
    ## add_header_above(c(' ' = 1, 'Raw' = 2, '20' = 2, '30' = 2)) %>%
    ## save_kable('../plots/controls/tables/vary_qual_one_plasmid_pct_non_canonical_bcs.html')

### - Combine Quality and Clustering Error Percentages --------------------

combo_one_plasmid_qual_cluster_counts = bind_rows(
    rename(one_plasmid_counts, correction = cluster),
    filter(vary_qual_one_plasmid_counts, quality != 'Raw') %>% rename(correction = quality))

combo_one_plasmid_qual_cluster_counts %>%
    group_by(library, correction) %>%
    summarize(non_canonical_pct =
                  sum((barcode != ten_control_barcodes[['BC7_1']]) * count) /
                  sum(count)) %>%
    ungroup() %>%
    mutate(non_canonical_pct = scales::percent(non_canonical_pct, accuracy = 0.01)) %>% 
    pivot_wider(names_from = correction,
                values_from = c(non_canonical_pct)) %>%
    left_join(tibble(library = c('BC1_102', 'BC1_103', 'BC1_104', 'BC1_105'),
                     library_name = c('1P-A', '1P-B', '1P-C', '1P-D'))  # This is to rename libraries
        ) %>%
    select(library_name, Raw, `20`, `30`, L1, L2, L3) %>%
    kbl(#caption = 'Percent Non-Plasmid Barcodes By Quality Cutoff And Clustering',
        col.names = c('Library', 'Raw', '20', '30', 'L1', 'L2', 'L3'),
        align = 'lrrrrrr') %>%
    column_spec(2:7, width_min = '0.8in') %>%
    ## collapse_rows(columns = 1, valign = 'top') %>%
    kable_styling(bootstrap_options='striped', full_width=FALSE, html_font = 'Arial') %>%
    add_header_above(c(' ' = 2, 'Quality Cutoff' = 2, 'Clustering' = 3)) %>%
    save_kable('../plots/controls/tables/combo_one_plasmid_qual_cluster_pct_non_canonical_bcs.html')


## -- Ten-Plasmid Controls ------------------------------------------------

vary_qual_ten_plasmid_counts = expand_grid(library = c('10diff1', '10diff2', '10diff3',
                                                       '1_01', '1_12', '10_10', '10_12'),
                                           quality = c('Raw', '20', '30')) %>%
    ## filter(cluster != 'raw') %>%
    mutate(data = map2(.x = library, .y = quality,
                       ~ if (.y != 'Raw') {
                             read_table(paste0('../barcodes_e0.1/cutadapt_fastqs/JA19375/',
                                               .x, '_qual_gte', .y, '_counts.txt'),
                                        col_names = c('count', 'barcode'),
                                        progress = FALSE)
                         }
                         else {
                             read_table(paste0('../barcodes_e0.1/JA19375/',
                                               .x, '_cutadapt_counts.txt'),
                                        col_names = c('count', 'barcode'),
                                        progress = FALSE)})) %>%
    unnest(., data) %>% 
    mutate(library = rename_10P[library])

vary_qual_selected_barcodes = vary_qual_ten_plasmid_counts %>%
    group_by(library, quality) %>%
    mutate(frac = count / sum(count)) %>%
    ungroup() %>%
    split(.$library) %>%
    imap(~ select(.x, - library)) %>%
    imap(~ filter(.x, ! barcode %in% ten_control_barcodes)) %>%
    imap(~ group_by(.x, quality)) %>%
    imap(~ top_n(.x, 2, wt = frac)) %>%
    imap(~ pull(.x, barcode)) %>%
    imap(~ union(.x, ten_control_barcodes)) %>%
    enframe(name = 'library', value = 'selected_barcodes') %>%
    unnest(cols = selected_barcodes)

vary_qual_ten_plasmid_counts %>%
    complete(library, quality, barcode, fill = list(count = 0)) %>%
    group_by(library, quality) %>%
    mutate(frac = count / sum(count)) %>%
    ungroup() %>%
    right_join(selected_barcodes,
               by = c('library' = 'library', 'barcode' = 'selected_barcodes')) %>% 
    ## right_join(selected_barcodes %>% filter(library %in% main_text_10P_controls),
    ##            by = c('library' = 'library', 'barcode' = 'selected_barcodes')) %>% 
    group_by(library) %>%
    nest() %>%
    mutate(plot = map(data, ~ ggplot(data = .x, aes(x = reorder(barcode, frac), y = frac)) +
                                geom_col(aes(fill = barcode)) +
                                coord_flip() +
                                facet_grid(cols = vars(factor(quality, levels = c('Raw', '20', '30')))) +
                                scale_fill_manual(name = 'Controls',
                                                  values = c(colors_ten_control_barcodes, 'gray85'),
                                                  limits = c(ten_control_barcodes, 'Other'),
                                                  na.value = 'gray85') +
                                scale_y_continuous('Percent Of Total Barcodes',
                                                   labels = scales::label_percent(accuracy = 1),
                                                   ## breaks = c(.125, 0.25, 0.375, 0.5),
                                                   minor_breaks = NULL) +
                                labs(tag = library) +
                                theme_minimal() +
                                theme(panel.spacing.y = unit(2, "lines"),
                                      axis.title = element_blank(),
                                      axis.text.y = element_text(family = 'IBM_Plex_Mono'),
                                      legend.text = element_text(family = 'IBM_Plex_Mono'),
                                      legend.position = 'none'
                                )
                      )
           ) %>%
    pull(plot) %>%
    reduce(`/`) +
    ## plot_annotation(
    ##     tag_levels = 'A',
    ##     tag_prefix = 'Fig. ',
    ##     theme = theme(plot.title = element_text(size = 18))
    ## ) +
    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->
    p_vary_qual_ten_plasmids
    ## xlab(label = 'Percent Of Total Barcodes') ->
    ## p_10P_ADG_L_params_pct_bar
    ## p_10P_A_D_G_L_params_pct_bar
    ## p_10P_C_D_L_params_pct_bar
    ## p_10P_L_params_pct_bar

ggsave('vary_qual_controls_all_10P_sqrt_bar.pdf',
       path = '../plots/controls',
       plot = p_vary_qual_ten_plasmids &
           scale_y_continuous(trans = 'sqrt',
                              labels = scales::label_percent(accuracy = 1),
                              breaks = c(0.01, 0.10, 0.25, 0.5, 1),
                              ## n.breaks = 5,
                              minor_breaks = NULL),
       width = 15,
       height = 30)

ggsave('vary_qual_controls_all_10P_bar.pdf',
       path = '../plots/controls',
       plot = p_vary_qual_ten_plasmids,
       width = 15,
       height = 30)

## -- Examine Read Counts For Individual 10P Controls ---------------------

vary_qual_ten_plasmid_counts %>%
    right_join(selected_barcodes, by = c('library' = 'library', 'barcode' = 'selected_barcodes')) %>%
    filter(library == '10P-D') %>%
    select(-library) %>%
    mutate(count = scales::comma(count, accuracy = 1)) %>%
    pivot_wider(names_from = quality, values_from = count) %>%
    kbl(caption = '10P-D: Num Barcodes By Quality Cutoff',
        col.names = c('Library', '# barcodes', '# barcodes', '# barcodes'),
        align = 'lrrr') %>%
    column_spec(1, bold = TRUE) %>%
    ## collapse_rows(columns = 1, valign = 'top') %>%
    kable_styling(bootstrap_options='striped', full_width=FALSE) %>%
    ## add_header_above(c(' ' = 1, 'Raw' = 2, '20' = 2, '30' = 2)) %>%
    save_kable('../plots/controls/tables/vary_qual_10P_num_barcodes.html')
