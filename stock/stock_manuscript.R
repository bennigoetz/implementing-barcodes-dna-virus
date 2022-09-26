library(tidyverse)
library(ComplexUpset)
library(patchwork)


# --- Select Colors -------------------------------------------------------

colors_ten_control_barcodes = colorspace::sequential_hcl(10, palette = 'Hawaii')
names(colors_ten_control_barcodes) = as.vector(ten_control_barcodes)

# --- Counts --------------------------------------------------------------

raw_dir_JA19161 = file.path('..', 'barcodes_e0.1', 'JA19161')
raw_dir_JA19375 = file.path('..', 'barcodes_e0.1', 'JA19375')

dir_JA19161 = file.path('..', 'barcodes_e0.1', 'message_passing_L3', 'JA19161')
dir_JA19375 = file.path('..', 'barcodes_e0.1', 'message_passing_L3', 'JA19375')

## -- Clustered Stock -----------------------------------------------------

stock_cluster_files = list.files(c(dir_JA19161, dir_JA19375),
                               full.names=TRUE,
                               pattern='BC7[_|b].*mp_L3_clusters.tsv')

stock_clustered_counts = stock_cluster_files %>%
    map(~read_tsv(., col_names=c('barcode',
                                 str_extract(., 'BC7[^_]*'),
                                 'elements'))) %>%
    map(~select(., -elements)) %>%
    reduce(full_join) %>%
    replace(is.na(.), 0) %>%
    pivot_longer(-barcode, names_to='BC7_run', values_to='count')

stock_clustered_counts %<>%
    group_by(barcode) %>%
    summarize(count = sum(count))

## -- Plasmid and Ligand --------------------------------------------------

plasmid_raw_counts = read_table(file.path(raw_dir_JA19375, 'BC7P_cutadapt_counts.txt'),
                                col_names = c('BC7P', 'barcode')) %>%
    mutate(barcode = replace_na(barcode, ''))

plasmid_clustered_counts = read_tsv(file.path(dir_JA19375, 'BC7P_mp_L3_clusters.tsv'),
                                    col_names = c('barcode', 'count', 'elements')) %>%
    select(- elements) 

## -- Raw Stock -----------------------------------------------------------

stock_raw_files = list.files(c(raw_dir_JA19161, raw_dir_JA19375),
                               full.names=TRUE,
                               pattern='BC7[_|b].*cutadapt_counts.txt')

stock_raw_counts = stock_raw_files %>%
    map(~read_table(., col_names=c(str_extract(., 'BC7[^_]*'),
                                   'barcode'))) %>%
    map(~mutate(., barcode = replace_na(barcode, ''))) %>%
    reduce(full_join) %>%
    replace(is.na(.), 0) %>%
    pivot_longer(-barcode, names_to='BC7_run', values_to='count')

stock_raw_counts %<>%
    group_by(barcode) %>%
    summarize(count = sum(count))

stock_raw_counts %>%
    filter(barcode != '') %>%
    write_tsv('BC7_combined_stock_counts.tsv',
              col_names = FALSE)

stock_L3_counts = read_tsv('BC7_combined_stock_mp_L3_clusters.tsv',
                           col_names = c('barcode', 'count', 'elements'),
                           progress = FALSE) %>%
    select(-elements)
                            

## -- Plot Count Rankings -------------------------------------------------

### - Stock ---------------------------------------------------------------

stock_barcode_ranking = stock_clustered_counts %>%
    arrange(desc(count)) %>%
    rownames_to_column('ranking') %>%
    mutate_at('ranking', as.integer)

stock_barcode_ranking %>%
    filter(ranking <= 15000) %>%
    ggplot(aes(ranking, count)) +
    geom_line(color = 'cornflowerblue') +
    scale_x_continuous(
        name='Rank'
        ) +
    scale_y_continuous(
        'Counts',
        labels=scales::label_number_si()
        ) +
    labs(title = 'Ranked Stock Counts',
         subtitle = 'Linear Scale') +
    theme_minimal() -> p_stock_count_ranking

ggsave('man_stock_count_ranking.pdf',
       path='../plots/stock_mp',
       plot = p_stock_count_ranking)
    
ggsave('man_stock_count_ranking_log10.pdf',
       path='../plots/stock_mp',
       plot = p_stock_count_ranking +
           scale_y_log10('Counts', labels=scales::label_number_si()) +
           labs(subtitle = 'Log 10 Scale'))

### - Raw Stock Counts ----------------------------------------------------

stock_raw_barcode_ranking = stock_raw_counts %>%
    arrange(desc(count)) %>%
    rownames_to_column('ranking') %>%
    mutate_at('ranking', as.integer)

stock_raw_barcode_ranking %>%
    filter(ranking <= 15000) %>%
    ggplot(aes(ranking, count)) +
    geom_line(color = 'cornflowerblue') +
    scale_x_continuous(
        name='Rank'
        ) +
    scale_y_continuous(
        'Counts',
        labels=scales::label_number_si()
        ) +
    labs(title = 'Ranked Raw Stock Counts',
         subtitle = 'Linear Scale') +
    theme_minimal() -> p_stock_raw_count_ranking

ggsave('man_stock_raw_count_ranking.pdf',
       path='../plots/stock_mp',
       plot = p_stock_raw_count_ranking)
    
ggsave('man_stock_raw_count_ranking_log10.pdf',
       path='../plots/stock_mp',
       plot = p_stock_raw_count_ranking +
           scale_y_log10('Counts', labels=scales::label_number_si()) +
           labs(subtitle = 'Log 10 Scale'))


### - Plasmid -------------------------------------------------------------

plasmid_barcode_ranking = plasmid_clustered_counts %>%
    arrange(desc(count)) %>%
    rownames_to_column('ranking') %>%
    mutate_at('ranking', as.integer)

plasmid_barcode_ranking %>%
    ## filter(ranking <= 15000) %>%
    ggplot(aes(ranking, count)) +
    geom_line(color = 'cornflowerblue') +
    scale_x_continuous(
        name='Rank'
        ) +
    scale_y_continuous(
        'Counts',
        labels=scales::label_number_si()
        ) +
    labs(title = 'Ranked Plasmid Counts',
         subtitle = 'Linear Scale') +
    theme_minimal() -> p_plasmid_count_ranking

ggsave('man_plasmid_count_ranking.pdf',
       path='../plots/stock_mp',
       plot = p_plasmid_count_ranking)
    
ggsave('man_plasmid_count_ranking_log10.pdf',
       path='../plots/stock_mp',
       plot = p_plasmid_count_ranking +
           scale_y_log10('Counts', labels=scales::label_number_si()) +
           labs(subtitle = 'Log 10 Scale'))

### - Stock vs Plasmid ----------------------------------------------------

bind_rows(plasmid_barcode_ranking %>%
            select(- barcode) %>%
            mutate(library = 'plasmid') %>%
            mutate(pct = count / sum(count)),
          stock_barcode_ranking %>%
            select(- barcode) %>%
            mutate(library = 'stock') %>%
            mutate(pct = count / sum(count))
          ) %>%
    filter(ranking <= 5965) %>%  # Max rank of plasmid ranking; an automated method here would be better
    ggplot(aes(x = ranking, y = pct, color = as_factor(library))) +
    geom_line() +
    labs(title = 'Stock vs Plasmid Ranked Percent (Clustered)') +
    theme_minimal() -> p_plasmid_stock_pct_ranking

ggsave('man_plasmid_stock_pct_ranking.pdf',
       path = '../plots/stock_mp')

ggsave('man_plasmid_stock_pct_ranking_log10.pdf',
       path = '../plots/stock_mp',
       plot = p_plasmid_stock_pct_ranking +
           scale_y_log10('Percent'))
       
### - Stock Raw vs Clustered ----------------------------------------------

## Percent
bind_rows(stock_barcode_ranking %>%
            select(-barcode) %>%
            mutate(library = 'clustered') %>%
            mutate(pct = count / sum(count)),
          stock_raw_barcode_ranking %>%
            select(-barcode) %>%
            mutate(library = 'raw') %>%
            mutate(pct = count /sum(count))
          ) %>%
    filter(ranking <= 15000) %>%
    ggplot(aes(x = ranking, y = pct, color = as_factor(library))) +
    geom_line() +
    labs(title = 'Stock Clustered vs Raw Percent') +
    theme_minimal() -> p_stock_raw_clustered_pct_ranking

ggsave('stock_raw_clustered_pct_ranking.pdf',
       path = '../plots/stock_mp')

ggsave('stock_raw_clustered_pct_ranking_log10.pdf',
       path = '../plots/stock_mp',
       plot = p_stock_raw_clustered_pct_ranking +
           scale_y_log10('Percent'))
       
## Counts
bind_rows(stock_barcode_ranking %>%
            select(-barcode) %>%
            mutate(library = 'clustered'),
          stock_raw_barcode_ranking %>%
            select(-barcode) %>%
            mutate(library = 'raw')
          ) %>%
    filter(ranking <= 15000) %>%
    ggplot(aes(x = ranking, y = count, color = as_factor(library))) +
    geom_line() +
    labs(title = 'Stock Clustered vs Raw Counts') +
    theme_minimal() -> p_stock_raw_clustered_count_ranking

ggsave('stock_raw_clustered_count_ranking.pdf',
       path = '../plots/stock_mp')

ggsave('stock_raw_clustered_count_ranking_log10.pdf',
       path = '../plots/stock_mp',
       plot = p_stock_raw_clustered_count_ranking +
           scale_y_log10('Percent'))
       

### - Plasmid Raw vs Clustered --------------------------------------------

plasmid_raw_barcode_ranking = plasmid_raw_counts %>%
    rename(count = BC7P) %>%
    arrange(desc(count)) %>%
    rownames_to_column('ranking') %>%
    mutate_at('ranking', as.integer)

## Counts
bind_rows(plasmid_barcode_ranking %>%
            select(-barcode) %>%
            mutate(library = 'clustered'),
          plasmid_raw_barcode_ranking %>%
            select(-barcode) %>%
            mutate(library = 'raw')
          ) %>%
    filter(ranking <= 10000) %>%
    ggplot(aes(x = ranking, y = count, color = as_factor(library))) +
    geom_line() +
    labs(title = 'Plasmid Clustered vs Raw Counts') +
    theme_minimal() -> p_plasmid_raw_clustered_count_ranking

ggsave('plasmid_raw_clustered_count_ranking.pdf',
       path = '../plots/stock_mp',
       plot = p_plasmid_raw_clustered_count_ranking)

ggsave('plasmid_raw_clustered_count_ranking_log10.pdf',
       path = '../plots/stock_mp',
       plot = p_plasmid_raw_clustered_count_ranking +
           scale_y_log10('Counts, Log 10'))
       

## Compare raw plasmid rankings with various distances for Levenshtein clustering

cluster_names = c('L1', 'L2', 'L3')
names(cluster_names) = cluster_names

clustering_type_colors = c(
    'Raw' = 'darkgoldenrod3',
    'L1' = 'lightskyblue1',
    'L2' = 'skyblue',
    'L3' = 'cornflowerblue'
    )

clustering_type_colors[2:4] %>%
    ## imap(~ file.path('..', 'barcodes_e0.1', paste('message_passing', .y, sep='_'), 'JA19375')) %>%
    imap(~ read_tsv(file.path(paste('BC7P_mp', .y, 'clusters.tsv', sep='_')),
                    col_names = c('barcode', 'count', 'elements'))) %>%
    imap(~ mutate(.x, library = .y)) %>%
    map(~ arrange(., desc(count))) %>%
    map(~ rownames_to_column(., 'ranking')) %>%
    map(~ mutate_at(., 'ranking', as.integer)) %>%
    map(~ select(., -elements, -barcode)) %>%
    reduce(bind_rows) %>%
    bind_rows(plasmid_raw_barcode_ranking %>%
              select(- barcode) %>%
              mutate(library = 'Raw')) %>%
    mutate(library = as_factor(library)) %>%
    filter(ranking <= 10000) %>%
    ggplot(aes(x = ranking, y = count, color = library)) +
    geom_line() +
    scale_color_manual(name = 'Clustering',
                       values = clustering_type_colors,
                       breaks = c('Raw', 'L1', 'L2', 'L3')) +
    labs(#title = 'Counts Of Ranked Barcodes By Clustering Parameters',
         x = 'Rank of Centroid By Cluster Size') +
    theme_minimal() +
    theme(
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 2),
        ) -> p_plasmid_raw_L1_L3_count_ranking

ggsave('plasmid_raw_L1_L3_count_ranking_log10.pdf',
       path = '../plots/stock',
       plot = p_plasmid_raw_L1_L3_count_ranking +
           scale_y_log10('Read Counts In Cluster, Log 10'))

## Compare raw viral stock rankings with various distances for Levenshtein clustering
## (same as plasmid thing above)

stock_raw_counts = read_tsv('BC7_combined_stock_counts.tsv',
                            col_names = c('barcode', 'count'))

stock_raw_barcode_ranking = stock_raw_counts %>%
    mutate(library = 'Raw') %>%
    mutate(ranking = row_number(desc(count)))

stock_L1_to_L3_cluster_counts =
    tibble(library = c('L1', 'L2', 'L3'),
           filenames = c('BC7_combined_stock_mp_L1_clusters.tsv',
                         'BC7_combined_stock_mp_L2_clusters.tsv',
                         'BC7_combined_stock_mp_L3_clusters.tsv')) %>%
    mutate(data = map(.x = filenames,
                       ~ read_tsv(.x,
                                col_names = c('barcode', 'count', 'elements'),
                                progress = FALSE) %>%
                       select(-elements)
                       )
           ) %>%
    unnest(., data) %>%
    select(- filenames) %>%
    group_by(library) %>%
    ## arrange(desc(count), .by_group = TRUE) %>%
    mutate(ranking = row_number(desc(count))) %>%
    ungroup()

stock_raw_barcode_ranking %>%
    bind_rows(stock_L1_to_L3_cluster_counts) %>%
    filter(ranking <= 10000) %>%
    ggplot(aes(x = ranking, y = count, color = library)) +
    geom_line() +
    scale_color_manual(name = 'Clustering',
                       values = clustering_type_colors,
                       breaks = c('Raw', 'L1', 'L2', 'L3')) +
    scale_y_log10('Counts, Log 10') +
    labs(#title = 'Counts Of Ranked Barcodes By Clustering Parameters',
        x = 'Ranking') +
    theme_minimal()

ggsave('stock_raw_to_L3_ranking.pdf',
       path = '../plots/stock')
    

# --- Another data frame for BC7P clusters --------------------------------

plasmid_cluster_counts =
    tibble(cluster = c('Raw', 'L1', 'L2', 'L3')) %>%
    mutate(data = map(.x = cluster,
                       ~ if (.x != 'Raw') {
                             read_tsv(paste('BC7P_mp', .x, 'clusters.tsv', sep = '_'),
                                      col_names = c('barcode', 'count', 'elements'),
                                      progress = FALSE) %>%
                                 select(-elements)
                         }
                         else {
                             read.table('../barcodes_e0.1/JA19375/BC7P_cutadapt_counts.txt',
                                        header = FALSE,
                                        col.names = c('count', 'barcode'),
                                        fill = TRUE,
                                        stringsAsFactors = FALSE)
                               })) %>%
    unnest(., data)

## -- Set Intersections, Upset Comparing Cluster Parameters ----------------
    
plasmid_cluster_counts %>%
    group_by(barcode) %>%
    mutate(total_count = sum(count)) %>%
    ungroup() %>%
    mutate(intersection = 1) %>%
    pivot_wider(id_cols = c(cluster, barcode, total_count), names_from = cluster, values_from = intersection) %>%
    replace_na(list(Raw = 0, L1 = 0, L2 = 0, L3 = 0)) %>%
    upset(c('Raw', 'L1', 'L2', 'L3'),
          name = 'Intersection Of Clustered Barcodes',
          set_sizes = FALSE,
          annotations = list(
              'Counts' = ggplot(mapping = aes(x = intersection, weight = total_count)) +
                  geom_bar(fill = 'cornflowerblue') +   # weight in aes makes geom_bar sum total_count
                  scale_y_continuous(name = 'Intersection count')
              )
          ) -> p_plasmid_clusters_upset

ggsave('plasmid_clusters_upset.pdf',
       path = '../plots/stock_mp',
       plot = p_plasmid_clusters_upset,
       width = 10)

# --- Comparing Plasmid, Ligated, Viral Libraries ---------------------------

L3_plas_lig_viral_cluster_counts =
    tibble(library = c('Plasmid', 'Ligated', 'Viral'),
           filenames = c('BC7P_mp_L3_clusters.tsv',
                         'BC7lig_mp_L3_clusters.tsv',
                         'BC7_combined_stock_mp_L3_clusters.tsv')) %>%
    mutate(data = map(.x = filenames,
                       ~ read_tsv(.x,
                                col_names = c('barcode', 'count', 'elements'),
                                progress = FALSE) %>%
                       select(-elements)
                       )
           ) %>%
    unnest(., data) %>%
    select(- filenames)

L3_cum99pct_plas_lig_viral_cluster_counts =
    L3_plas_lig_viral_cluster_counts %>%
    group_by(library) %>%
    mutate(total = sum(count),
           pct = count / total,
           cum_pct = cumsum(pct)) %>%
    filter(cum_pct <= 0.99) %>%
    ungroup()

### Barcodes shorter than 12bp

L3_plas_lig_viral_cluster_counts %>%
    mutate(length = str_length(barcode)) %>%
    group_by(barcode) %>%
    mutate(total_count = sum(count)) %>%
    ungroup() %>%
    group_by(library) %>%
    mutate(total = sum(count),
           pct = count /total,
           cum_pct = cumsum(pct)) %>%
    filter(cum_pct <= 0.99) %>%
    ungroup() %>%
    filter(length < 12) %>%
    mutate(library = if_else(library == 'Viral', 'Virus', library)) %>%  # Necessary to eplace 'Viral' with 'Virus' in plot
    mutate(library = if_else(library == 'Ligated', 'Ligated virus genomes', paste(library, 'library'))) %>%
    mutate(intersection = 1) %>%
    pivot_wider(id_cols = c(library, barcode, total_count, length),
                names_from = library,
                values_from = intersection) %>%
    replace_na(list(Plasmid = 0, Ligated = 0, Viral = 0)) %>%
    upset(c('Virus library', 'Ligated virus genomes', 'Plasmid library'),
          name = 'Intersection of Barcodes Shorter Than 12nt Between Libraries',
          set_sizes = FALSE,
          sort_sets = FALSE,
          sort_intersections = 'ascending',
          sort_intersections_by = 'degree',
          intersections = list(
              'Virus library',
              c('Virus library', 'Plasmid library'),
              c('Virus library', 'Ligated virus genomes'),
              c('Virus library', 'Ligated virus genomes', 'Plasmid library')
              ),
          annotations = list(
              'Counts' = ggplot(mapping = aes(x = intersection, weight = total_count)) +
                  geom_bar(fill = 'cornflowerblue') +   # weight in aes makes geom_bar sum total_count
                  scale_y_continuous(name = 'Intersection count')
              )
          ) -> p_L3_plas_lig_viral_lt12bp_upset

## ggsave('L3_plas_lig_viral_lt12bp_upset.pdf',
##        path = '../plots/stock_mp',
##        plot = p_L3_plas_lig_viral_lt12bp_upset,
##        width = 10)

ggsave('L3_plas_lig_viral_lt12bp_cum99pct_upset.pdf',
       path = '../plots/stock',
       plot = p_L3_plas_lig_viral_lt12bp_upset,
       width = 10)

### Viral counts found in plasmid and ligated as well

L3_plas_lig_viral_cluster_counts %>%
    # 99% cum count cutoff
    group_by(library) %>%
    mutate(total = sum(count),
           pct = count / total,
           cum_pct = cumsum(pct)) %>%
    filter(cum_pct <= 0.99) %>%
    ungroup() %>%
    # Create column where the same barcodes all have viral count, to prepare for pivot
    mutate(count = if_else(library == 'Viral', count, 0)) %>%
    group_by(barcode) %>%
    mutate(viral_count = sum(count)) %>%
    ungroup() %>%
    mutate(library = if_else(library == 'Viral', 'Virus', library)) %>%
    mutate(library = if_else(library == 'Ligated', 'Ligated virus genomes', paste(library, 'library'))) %>%
    mutate(intersection = 1) %>%
    pivot_wider(id_cols = c(library, barcode, viral_count),
                names_from = library,
                values_from = intersection) %>%
    replace_na(list(`Plasmid library`= 0, `Ligated virus genomes` = 0, `Viral library` = 0)) %>%
    upset(c('Virus library', 'Ligated virus genomes', 'Plasmid library'),
          name = 'Intersection of Barcodes Between Libraries',
          set_sizes = FALSE,
          sort_sets = FALSE,
          sort_intersections = 'ascending',
          sort_intersections_by = 'degree',
          intersections = list(
              c('Plasmid library', 'Virus library'),
              c('Ligated virus genomes', 'Virus library'),
              c('Plasmid library', 'Ligated virus genomes', 'Virus library'),
              'Virus library'
              ),
          ## intersections = list(
          ##     'Viral',
          ##     c('Viral', 'Plasmid'),
          ##     c('Viral', 'Ligated'),
          ##     c('Viral', 'Ligated', 'Plasmid')
          ##     ),
          annotations = list(
              'Counts' = ggplot(mapping = aes(x = intersection, weight = viral_count)) +
                  geom_bar(fill = 'cornflowerblue') +   # weight in aes makes geom_bar sum total_count
                  scale_y_continuous(name = 'Intersection count')
              )
          ) -> p_L3_viral_counts_intersect_plas_lig_upset

ggsave('L3_cum99pct_viral_counts_in_plas_lig_upset.pdf',
       path = '../plots/stock',
       plot = p_L3_viral_counts_intersect_plas_lig_upset, 
       width = 10)

## Combine the viral, ligated, plasmid UpSet plots into a patchwork

ggsave('L3_cum99pct_plas_lig_viral_all_lt12bp_combo_upset.pdf',
       path = '../plots/stock',
       plot = (wrap_elements(p_L3_viral_counts_intersect_plas_lig_upset) /
               wrap_elements(p_L3_plas_lig_viral_lt12bp_upset)) +
           plot_layout(nrow = 2) +
           plot_annotation(tag_levels = 'A'),
       width = 10,
       height = 14)

ggsave('L3_cum99pct_plas_lig_viral_all_lt12bp_combo_upset.pdf',
       path = '../plots/stock',
       plot = wrap_plots(p_L3_viral_counts_intersect_plas_lig_upset,
                         p_L3_plas_lig_viral_lt12bp_upset,
                         ncol = 1) +
      plot_annotation(tag_levels = c('A', '')),
           ## plot_annotation(tag_levels = 'A') +
           ## plot_layout(nrow = 6),
       width = 10,
       height = 14)


## Calculate counts, this replicates the code producing the intersection table.

L3_cum99pct_viral_intersection_table = L3_cum99pct_plas_lig_viral_cluster_counts %>%
    mutate(count = if_else(library == 'Viral', count, 0)) %>%
    group_by(barcode) %>%
    mutate(viral_count = sum(count)) %>%
    ungroup() %>%
    mutate(intersection = 1) %>%
    pivot_wider(id_cols = c(library, barcode, viral_count), names_from = library, values_from = intersection) %>%
    replace_na(list(`Plasmid`= 0, `Ligated` = 0, `Viral` = 0))

## Calculate percentage of counts in plasmid, ligated, viral intersection
(L3_cum99pct_viral_intersection_table %>% filter(Plasmid == 1, Ligated == 1, Viral == 1) %>% pull(viral_count) %>% sum()) / (L3_cum99pct_viral_intersection_table %>% filter(Viral == 1) %>% pull(viral_count) %>% sum())

## Calculate percentage of counts in just viral 
(L3_cum99pct_viral_intersection_table %>% filter(Plasmid == 0, Ligated == 0, Viral == 1) %>% pull(viral_count) %>% sum()) / (L3_cum99pct_viral_intersection_table %>% filter(Viral == 1) %>% pull(viral_count) %>% sum())

## Same thing, but for <12 bp
L3_cum99pct_12bp_viral_intersection_table = L3_cum99pct_plas_lig_viral_cluster_counts %>%
    mutate(length = str_length(barcode)) %>%
    filter(length < 12) %>%
    mutate(count = if_else(library == 'Viral', count, 0)) %>%
    group_by(barcode) %>%
    mutate(viral_count = sum(count)) %>%
    ungroup() %>%
    mutate(intersection = 1) %>%
    pivot_wider(id_cols = c(library, barcode, viral_count), names_from = library, values_from = intersection) %>%
    replace_na(list(`Plasmid`= 0, `Ligated` = 0, `Viral` = 0))

## Calculate percentage of counts in plasmid, ligated, viral intersection
(L3_cum99pct_12bp_viral_intersection_table %>% filter(Plasmid == 1, Ligated == 1, Viral == 1) %>% pull(viral_count) %>% sum()) / (L3_cum99pct_12bp_viral_intersection_table %>% filter(Viral == 1) %>% pull(viral_count) %>% sum())

## Calculate percentage of counts in just viral 
(L3_cum99pct_12bp_viral_intersection_table %>% filter(Plasmid == 0, Ligated == 0, Viral == 1) %>% pull(viral_count) %>% sum()) / (L3_cum99pct_12bp_viral_intersection_table %>% filter(Viral == 1) %>% pull(viral_count) %>% sum())

### Correlation between plasmid and viral abundance

L3_cum99pct_viral_and_plasmid_bcs = intersect(
    L3_cum99pct_plas_lig_viral_cluster_counts %>%
    filter(library == 'Viral') %>%
    pull(barcode),
    L3_cum99pct_plas_lig_viral_cluster_counts %>%
    filter(library == 'Plasmid') %>%
    pull(barcode)
    )

L3_cum99pct_plas_lig_viral_cluster_counts %>%
    filter(library != 'Ligated' & barcode %in% L3_cum99pct_viral_and_plasmid_bcs) %>%
    mutate(pct = 100 * pct) %>%
    pivot_wider(id_cols = c(library, barcode, pct),
                names_from = library,
                values_from = pct) %>%
    ggplot(aes(x = Plasmid, y = Viral)) +
    geom_point(color = 'cornflowerblue', alpha = 0.5) +
    ## geom_smooth(method = 'lm', se = 0, color='gray85') +
    geom_abline(color = 'gray85') +
    ## ylim(0, 0.006) +
    scale_y_log10(name = 'Log 10 Percent of Total Barcode Count In Virus Library') +
    scale_x_log10(name = 'Log 10 Percent of Total Barcode Count In Plasmid Library') +
    ## scale_y_log10(name = 'Percent Of Virus Library', labels = NULL) +
    ## scale_x_log10(name = 'Percent Of Plasmid Library', labels = NULL) +
    ## scale_y_log10(limits = c(0.002, 1)) +
    ## scale_x_log10(limits = c(0.002, 1)) +
    theme_minimal() +
    theme(
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 1)
        ) ->
    p_L3_cum99pct_viral_plasmid_pct_scatter_w_abline

ggsave('L3_cum99pct_viral_plasmid_pct_scatter.pdf',
       path = '../plots/stock_mp') 
       ## plot = p_L3_viral_counts_intersect_plas_lig_upset, 
       ## width = 10)

ggsave('L3_cum99pct_viral_plasmid_pct_w_abline_scatter.pdf',
       path = '../plots/stock', 
       plot = p_L3_cum99pct_viral_plasmid_pct_scatter_w_abline)
       ## width = 10)

## Calculate R^2 for regression, and correlation between variables
L3_cum99pct_plas_v_lig_pct_wide = L3_cum99pct_plas_lig_viral_cluster_counts %>%
    filter(library != 'Ligated' & barcode %in% L3_cum99pct_viral_and_plasmid_bcs) %>%
    mutate(pct = 100*pct) %>%
    pivot_wider(id_cols = c(library, barcode, pct), names_from = library, values_from = pct)

lm(Viral ~ Plasmid, data = L3_cum99pct_plas_v_lig_pct_wide) %>% summary()

# Spearman
cor(L3_cum99pct_plas_v_lig_pct_wide$Viral, L3_cum99pct_plas_v_lig_pct_wide$Plasmid, method = 'spearman')

# Pearson
cor(L3_cum99pct_plas_v_lig_pct_wide$Viral, L3_cum99pct_plas_v_lig_pct_wide$Plasmid, method = 'pearson')


## -- Examining Count Cutoffs ---------------------------------------------

L3_plas_lig_viral_cluster_counts %>%
    group_by(library) %>%
    arrange(desc(count))

## -- Constituent Barcodes Of Clusters-------------------------------------

stock_clusters = stock_cluster_files %>%
    map(~read_tsv(., col_names=c('barcode',
                                 str_extract(., 'BC7[^_]*'),
                                 'elements'))) %>%
    reduce(full_join) %>%
    replace(is.na(.), 0) %>%
    pivot_longer(-barcode, names_to='BC7_run', values_to='count')

stock_clusters %<>%
    group_by(barcode) %>%
    summarize(count = sum(count))

### - Do This Just For BC7 For Simplicity (and because I screwed the pooch on separate clustering)

BC7_clusters = read_tsv("../barcodes_e0.1/message_passing_L3/JA19161/BC7_mp_L3_clusters.tsv",
                        col_names=c('centroid', 'cluster_count', 'barcode')) %>%
    replace(is.na(.), 0)

BC7_raw_counts = read_table(file.path(raw_dir_JA19161, 'BC7_cutadapt_counts.txt'),
                                col_names = c('BC7', 'barcode')) %>%
    mutate(barcode = replace_na(barcode, ''))

BC7_clusters %>%
    separate_rows(barcode) %>%
    left_join(BC7_raw_counts) %>%
    filter(barcode == centroid) %>%
    arrange(desc(cluster_count)) %>%
    rownames_to_column('ranking') %>%
    mutate_at('ranking', as.integer) %>%
    ggplot(aes(x = ranking, y = cluster_count)) +
    geom_area(fill = 'grey90') +
    geom_area(aes(y = BC7), fill = 'cornflowerblue') +
    labs(title = 'Centroid Counts (Blue) Over Cluster Counts(Gray)',
         subtitle = 'Linear Scale') +
    theme_minimal() -> p_BC7_centroid_over_cluster_counts

ggsave('BC7_centroid_over_cluster_count_ranked.pdf',
       path = '../plots/stock_mp',
       plot = p_BC7_centroid_over_cluster_counts)

ggsave('BC7_centroid_over_cluster_count_ranked_log10.pdf',
       path = '../plots/stock_mp',
       plot = p_BC7_centroid_over_cluster_counts +
           labs(subtitle = 'Log 10 Scale') +
           scale_y_log10())

BC7_clusters %>%
    separate_rows(barcode) %>%
    left_join(BC7_raw_counts) %>%
    filter(barcode == centroid) %>%
    arrange(desc(cluster_count)) %>%
    rownames_to_column('ranking') %>%
    mutate_at('ranking', as.integer) %>%
    mutate(pct = BC7 / cluster_count) %>%
    ggplot(aes(x = ranking, y = pct)) +
    geom_point(color = 'gray90') +
    geom_smooth() +
    labs(title = 'Centroid Counts As Percentage Of Cluster Count') +
    theme_minimal() -> p_BC7_centroid_over_cluster_pct

ggsave('BC7_centroid_over_cluster_pct_ranked.pdf',
       path = '../plots/stock_mp',
       plot = p_BC7_centroid_over_cluster_pct)

## -- How Do Centroid Barcodes Overlap With Raw Barcodes? -----------------

(plasmid_barcode_ranking %>% pull(count) %>% sum()) * 0.99

plasmid_raw_barcode_ranking %>% filter(barcode %in% (plasmid_barcode_ranking %>% mutate(cum_count = cumsum(count)) %>% filter(cum_count < 10640497) %>% pull(barcode))) %>% pull(ranking) %>% max()

plasmid_raw_barcode_ranking %>% filter(ranking <= 5729) %>% filter(! barcode %in% plasmid_barcode_ranking$barcode) %>% pull(count) %>% sum()

plasmid_raw_barcode_ranking %>% filter(ranking <= 5729) %>%  pull(count) %>% sum()

plasmid_raw_barcode_ranking %>%
    filter(ranking <= 5729) %>%
    filter(! barcode %in% plasmid_barcode_ranking$barcode) %>%
    ggplot(aes(x = ranking, y = count)) +
    geom_area() +
    theme_minimal()

plasmid_raw_barcode_ranking %>%
    filter(ranking <= 10000) %>%
    mutate(library = if_else(! barcode %in% plasmid_barcode_ranking$barcode, 'grey70', 'cornflowerblue')) %>%
    ggplot(aes(x = ranking, y = count, fill = library)) +
    geom_area() +
    theme_minimal() -> p_plasmid_raw_bcs_among_top_L3_ranked

ggsave('plasmid_raw_bcs_among_top_L3_ranked.pdf',
       path = '../plots/stock_mp',
       plot = p_plasmid_raw_bcs_among_top_L3_ranked +
       scale_y_log10())

## Bin raw plasmid barcodes in groups of 500, calculate what percent are also centroids of L3 clustering

plasmid_raw_barcode_ranking %>%
    filter(ranking <= 10000) %>%
    mutate(bin = cut(ranking, seq(1, 10001, 500), right = FALSE, labels = FALSE)) %>%
    mutate(centroid = barcode %in% plasmid_barcode_ranking$barcode) %>%
    group_by(bin) %>%
    summarize(pct_centroid = sum(centroid) / 500) %>%
    ggplot(aes(x = bin, y = pct_centroid)) +
    geom_bar(stat = 'identity', fill = 'cornflowerblue') +
    scale_y_continuous(labels = scales::label_percent()) +
    labs(title = 'Plasmid: Similarity Between Raw Barcodes and L3 Centroids',
         subtitle = 'Percentage Of Raw Barcodes Also Centroids, Binned Into Groups Of 500',
         x = 'Bin',
         y = 'Percentage Of Raw Barcodes Also Centroids') +
    theme_minimal() -> p_plasmid_raw_barcodes_pct_centroid_bin500

ggsave('plasmid_raw_barcodes_also_centroid_bin500_bar.pdf',
       path = '../plots/stock_mp',
       plot = p_plasmid_raw_barcodes_pct_centroid_bin500,
       width = 10
       )

## -- What Percentage Of Clusters Do Centroids Make Up? -------------------

BC7_raw_counts = read_table(file.path(raw_dir_JA19161, 'BC7_cutadapt_counts.txt'),
                                col_names = c('BC7', 'barcode')) %>%
    mutate(barcode = replace_na(barcode, ''))

plasmid_L3_clusters = read_tsv('BC7P_mp_L3_clusters.tsv',
                               col_names=c('centroid', 'cluster_count', 'barcode')) %>%
    replace(is.na(.), 0) %>%
    separate_rows(barcode)


plasmid_cluster_counts %>%
    filter(cluster == 'Raw') %>%
    left_join(plasmid_L3_clusters) %>%
    filter(barcode == centroid) %>%
    arrange(desc(cluster_count)) %>%
    rownames_to_column('ranking') %>%
    mutate_at('ranking', as.integer) %>%
    mutate(pct = count / cluster_count) %>%
    ggplot(aes(x = ranking, y = pct)) +
    geom_point(color = 'gray85') +
    geom_smooth() +
    labs(x = 'Rank of Centroid By Cluster Size',
         y = 'Centroid Percentage of Cluster') +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::label_percent()) +
    theme_minimal() -> p_plasmid_centroid_pct_of_cluster_pct

ggsave('plasmid_centroid_pct_of_cluster_ranked.pdf',
       path = '../plots/stock',
       plot = p_plasmid_centroid_pct_of_cluster_pct)

## -- Testing Cumulative Percentage Cutoffs -------------------------------


tibble(cum_cutoff = c(0.9, 0.99, 0.999, 1)) %>%
    mutate(data = map(.x = cum_cutoff,
                      ~ (plasmid_cluster_counts %>%
                               filter(cluster == 'L3') %>%
                               mutate(total = sum(count),
                                      pct = count / total,
                                      cum_pct = cumsum(pct)) %>%
                         filter(cum_pct <= .x) %>%
                         mutate(bc_length = str_length(barcode)) %>%
                         group_by(bc_length) %>%
                         summarize(count = n_distinct(barcode))
                      )
                      )
           ) %>%
    unnest(., data) %>%
    mutate(cum_cutoff = scales::percent(cum_cutoff, accuracy = 0.1)) %>%
    mutate(cum_cutoff = factor(cum_cutoff, levels = c('90.0%', '99.0%', '99.9%', '100.0%'))) %>%
    ggplot(aes(x = bc_length, y = count)) +
    geom_bar(stat = 'identity',
             fill = 'cornflowerblue') +
    facet_wrap(vars(cum_cutoff), nrow = 2) +
    scale_x_continuous(
        'Barcode Length',
        breaks=seq(0, 75, 2),
        minor_breaks = NULL) +
    scale_y_sqrt(
        'Number of Barcodes',
        labels=scales::label_number_si(),
        minor_breaks = NULL) +
    theme_minimal() +
    theme(strip.background = element_rect(color = NA,
                                          fill="gray85"),
          panel.spacing = unit(7, "mm")) -> p_plasmid_L3_length_distribution_by_cum_cutoff

ggsave('plasmid_L3_cum_pct_various_length_distribution.pdf',
       path='../plots/stock_mp',
       plot = p_plasmid_L3_length_distribution_by_cum_cutoff,
       width = 10)


## -- Plasmids: Distance Between Centroids --------------------------------

plasmid_cluster_counts %>%
    filter(cluster == 'L3') %>%
    ## filter(count >= 100) %>%
    mutate(total = sum(count),
           pct = count /total,
           cum_pct = cumsum(pct)) %>%
    filter(cum_pct <= 0.99) %>%
    pull(barcode) %>%
    set_names(., nm=.) %>%
    adist() %>%
    as.data.frame() %>%
    rownames_to_column('stock_1') %>%
    as_tibble %>%
    pivot_longer(-stock_1, names_to='stock_2', values_to='dist') %>%
    filter(stock_1 < stock_2) %>%
    ggplot(., aes(dist)) +
    geom_histogram(binwidth = 1, fill = 'cornflowerblue') +
    labs(x = 'Levenshtein Distance',
         y = 'Count of Pairwise Distances',
         title = 'Plasmid Library') +
    scale_x_continuous(breaks = seq(0, 12, 4)) +
    scale_y_continuous(breaks = NULL) +
    theme_minimal() +
    theme(
        axis.text = element_text(size = 10)
        ) -> p_plasmid_L3_distance_distribution

ggsave('plasmid_L3_cum_pct_99_distance_distribution.pdf',
       path='../plots/stock',
       plot = p_plasmid_L3_distance_distribution)

plasmid_cluster_counts %>%
    filter(cluster == 'L3') %>%
    ## filter(count >= 100) %>%
    mutate(total = sum(count),
           pct = count /total,
           cum_pct = cumsum(pct)) %>%
    filter(cum_pct <= 0.99) %>%
    mutate(bc_length = str_length(barcode)) %>%
    group_by(bc_length) %>%
    summarize(count = n_distinct(barcode)) %>%
    ggplot(aes(x = bc_length, y = count)) +
    ## geom_bar(stat = 'identity',
    ##          fill = 'cornflowerblue') +
    geom_bar(stat = 'identity',
             fill = 'cornflowerblue') +
    scale_x_continuous(
        'Barcode length',
        breaks=seq(0, 75, 2),
        minor_breaks = NULL) +
    scale_y_continuous(
        'Count',
        labels=scales::label_number_si(),
        minor_breaks = NULL) +
    theme_minimal() -> p_plasmid_L3_length_distribution

ggsave('plasmid_L3_cum_pct_99_length_distribution.pdf',
       path='../plots/stock_mp',
       plot = p_plasmid_L3_length_distribution)

p_plasmid_L3_distance_length_distributions = 
    p_plasmid_L3_distance_distribution +
    p_plasmid_L3_length_distribution +
    plot_annotation(tag_levels = 'A')

## p_plasmid_L3_distance_length_distributions =
##     p_plasmid_L3_distance_length_distributions +
##     plot_annotation(tag_levels = 'A')

ggsave('plasmid_L3_cum_pct_99_distance-length_dist_bar.pdf',
       path = '../plots/stock_mp',
       plot = p_plasmid_L3_distance_length_distributions,
       height = 6,
       width = 14)

# --- UpSet Plots ---------------------------------------------------------

## -- Raw, Unclustered Barcodes -------------------------------------------

raw_dir_JA19161 = file.path('..', 'barcodes_e0.1', 'JA19161')
raw_dir_JA19375 = file.path('..', 'barcodes_e0.1', 'JA19375')

raw_stock_count_files = list.files(c(raw_dir_JA19161, raw_dir_JA19375),
                                   full.names=TRUE,
                                   pattern='BC7[_|b].*_counts.txt')

raw_stock_barcodes_separate = raw_stock_count_files %>%
    map(~ read_table(., col_names=c(str_extract(., 'BC7[^_]*'), 'barcode'))) %>%
    map(~ mutate(., barcode = replace_na(barcode, ''))) %>%
    set_names(c('Stock 1', 'Stock 2', 'Stock 3', 'Stock 4'))
    
raw_stock_barcodes = raw_stock_barcodes_separate %>%
    map(~ pull(., barcode)) %>%
    reduce(union) %>%
    list()
names(raw_stock_barcodes) = 'BC7'

## -- Compare L3 Viral Stock Runs -----------------------------------------

L3_stock_runs =
    tibble(library = c('Stock 1', 'Stock 2', 'Stock 3', 'Stock 4'),
           filenames = c('BC7_mp_L3_clusters.tsv',
                         'BC7bis1_mp_L3_clusters.tsv',
                         'BC7bis2_mp_L3_clusters.tsv',
                         'BC7bis3_mp_L3_clusters.tsv')) %>%
    mutate(data = map(.x = filenames,
                       ~ read_tsv(.x,
                                col_names = c('barcode', 'count', 'elements'),
                                progress = FALSE) %>%
                       select(-elements)
                       )
           ) %>%
    unnest(., data) %>%
    select(- filenames)

L3_stock_runs %>%
    group_by(library) %>%
    mutate(total = sum(count),
           pct = count / total,
           cum_pct = cumsum(pct)) %>%
    filter(cum_pct <= 0.99) %>%
    ungroup() %>%
    group_by(barcode) %>%
    mutate(total_count = sum(count)) %>%
    ungroup() %>%
    select(- count) %>%
    mutate(intersection = 1) %>%
    pivot_wider(id_cols = c(library, barcode, total_count),
                names_from = library,
                values_from = intersection) %>%
    replace_na(list(`Stock 1` = 0, `Stock 2` = 0, `Stock 3` = 0, `Stock 4` = 0)) %>%
    upset(c('Stock 4', 'Stock 3', 'Stock 2', 'Stock 1'),
          name = 'Intersection of Virus Library Technical Replicates',
          set_sizes = FALSE,
          sort_intersections = 'ascending',
          sort_intersections_by = 'degree',
          ## intersections = list(
          ##     'Viral',
          ##     c('Viral', 'Plasmid'),
          ##     c('Viral', 'Ligated'),
          ##     c('Viral', 'Ligated', 'Plasmid')
          ##     ),
          sort_sets = FALSE,
          annotations = list(
              'Counts' = ggplot(mapping = aes(x = intersection, weight = total_count)) +
                  geom_bar(fill = 'cornflowerblue') +   # weight in aes makes geom_bar sum total_count
                  scale_y_continuous(name = 'Intersection count')
              )
          ) -> p_stock_runs_L3_cutoff_upset

ggsave('stock_runs_L3_cum99pct_upset.pdf',
       path = '../plots/stock',
       plot = p_stock_runs_L3_cutoff_upset,
       width = 10)


# --- Calculate Distance Between Stock Barcodes Out Of Curiosity ----------

stock_barcodes = L3_cum99pct_plas_lig_viral_cluster_counts %>%
    filter(library == 'Viral') %>%
    pull(barcode) %>%
    set_names(., nm = .)

distance_between_stock_barcodes = adist(stock_barcodes) %>%
    as.data.frame() %>%
    rownames_to_column('barcode1') %>%
    as_tibble() %>%
    pivot_longer(-barcode1, names_to = 'barcode2', values_to = 'distance')

sphere_L3_stock_barcodes = read_tsv('BC7_sphere_combined_stock_L3_clusters.tsv',
                                    col_names = c('barcode', 'count', 'elements'),
                                    progress = FALSE) %>%
    select(-elements) %>%
    mutate(total = sum(count),
           pct = count / total,
           cum_pct = cumsum(pct)) %>%
    filter(cum_pct <= 0.99) %>%
    pull(barcode) %>%
    set_names(., nm = .)

distance_between_sphere_L3_stock_barcodes = adist(sphere_L3_stock_barcodes) %>%
    as.data.frame() %>%
    rownames_to_column('barcode1') %>%
    as_tibble() %>%
    pivot_longer(-barcode1, names_to = 'barcode2', values_to = 'distance')

## -- Distance To Nearest Stock Barcode -----------------------------------

distance_between_stock_barcodes %>%
    filter(barcode1 < barcode2) %>%
    group_by(barcode1) %>%
    summarize(min_dist = min(distance)) %>%
    arrange(min_dist) %>%
    count(min_dist) 

## L3 clustered
stock_L3_cum99pct_min_dists = L3_cum99pct_plas_lig_viral_cluster_counts %>%
    filter(library == 'Viral') %>%
    pull(barcode) %>%
    set_names(., nm = .) %>%
    adist(.) %>%
    as.data.frame() %>%
    rownames_to_column('barcode1') %>%
    as_tibble() %>%
    pivot_longer(- barcode1, names_to = 'barcode2', values_to = 'distance') %>%
    filter(barcode1 < barcode2) %>%
    group_by(barcode1) %>%
    summarize(min_dist = min(distance)) %>%
    count(min_dist) %>%
    arrange(min_dist)

stock_raw_rank_20k_mid_dists = stock_raw_barcode_ranking %>%
    filter(ranking <= 20000) %>%    # Can't allocate a vector big enough for all barcodes
    pull(barcode) %>%
    set_names(., nm = .) %>%
    adist(.) %>%
    as.data.frame() %>%
    rownames_to_column('barcode1') %>%
    as_tibble() %>%
    pivot_longer(- barcode1, names_to = 'barcode2', values_to = 'distance') %>%
    filter(barcode1 < barcode2) %>%
    group_by(barcode1) %>%
    summarize(min_dist = min(distance)) %>%
    count(min_dist) %>%
    arrange(min_dist)

stock_combined_raw_L3_cum99pct_min_dists =
    bind_rows(mutate(stock_L3_cum99pct_min_dists, barcodes = 'L3_cum99pct', pct = n / sum(n)),
              mutate(stock_raw_rank_20k_mid_dists, barcodes = 'Raw_20k', pct = n / sum(n))
        )

ggplot(stock_combined_raw_L3_cum99pct_min_dists,
       aes(x = min_dist)) +
    geom_col(aes(y = n), fill = 'cornflowerblue', position = 'identity') +
    facet_wrap(vars(barcodes)) +
    theme_minimal() -> p_min_dists_raw_L3_cum99pct

ggsave('min_dists_raw_20k_L3_cum99pct_bar.pdf',
       path = '../plots/stock',
       plot = p_min_dists_raw_L3_cum99pct)

ggplot(stock_combined_raw_L3_cum99pct_min_dists,
       aes(x = min_dist)) +
    geom_col(aes(y = pct), fill = 'cornflowerblue', position = 'identity') +
    facet_wrap(vars(barcodes)) +
    theme_minimal() -> p_min_dists_raw_L3_cum99pct_pct

ggsave('min_dists_raw_20k_L3_cum99pct_pct_bar.pdf',
       path = '../plots/stock',
       plot = p_min_dists_raw_L3_cum99pct_pct)


# --- Simulate Pairwise Distance Of Random 12mers -------------------------

set_of_distances_within_3 = c()

## for (trial in 1:20) {

random_bcs = c()
for (i in 1:5000) {
    twelve_mer = sample(c('A', 'C', 'G', 'T'), 12, replace = TRUE) %>% paste0(collapse="");
    random_bcs = c(random_bcs, twelve_mer)
}
random_bcs = set_names(random_bcs, nm = random_bcs)

sample_bc_vector <- function(bc_vector, num_samples, num_iterations, prob = NULL) {
    print(num_iterations)
    length_of_bc_samples = vector('integer', num_iterations)
    for (i in 1:num_iterations) {
        print(i)
        sampled_bcs = sample(bc_vector, num_samples, replace = FALSE, prob)
        print(sampled_bcs)
        length_of_bc_samples[i] = length(sampled_bcs)
        return(length_of_bc_samples)
        }
}

random_12mer_pairwise = adist(random_bcs) %>%
    as.data.frame() %>%
    rownames_to_column('barcode1') %>%
    as_tibble() %>%
    pivot_longer(-barcode1, names_to = 'barcode2', values_to = 'distance')
    ## Maybe need barcode1 < barcode2 or something so there's no double-counting pairs

distances_within_3 = random_12mer_pairwise %>%
    filter(distance > 0 & distance <= 3) %>%
    pull(distance) %>%
    length()

set_of_distances_within_3 = c(set_of_distances_within_3, distances_within_3)


average_of_pct_distances_within_3 = mean(set_of_distances_within_3) / 25000000

random_12mer_pairwise %>%
    filter(barcode1 < barcode2) %>%
    ggplot(., aes(distance)) +
    geom_histogram(binwidth = 1, fill = 'darkgoldenrod3') +
    labs(x = 'Levenshtein Distance',
         y = 'Count of Pairwise Distances',
         title = '5000 Random 12mers') +
    scale_x_continuous(breaks = seq(0, 12, 4)) +
    scale_y_continuous(breaks = NULL) +
    theme_minimal() +
    theme(
        axis.text = element_text(size = 10)
    ) -> p_random5k_distance_distribution

ggsave('random5k_distance_distribution.pdf',
       path='../plots/stock',
       plot = p_random5k_distance_distribution)

# Combine random distance distribution with the plasmid distribution with L3/cutoff

p_plasmid_L3_distance_and_random5k_distance_distribution =
    p_plasmid_L3_distance_distribution +
    p_random5k_distance_distribution +
    plot_annotation(tag_levels = 'A') &
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
          )

ggsave('plasmid_and_5krandom_distance_distribution.pdf',
       path = '../plots/stock',
       plot = p_plasmid_L3_distance_and_random5k_distance_distribution,
       width = 14)


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

    theme(plot.subtitle = element_blank(),
          plot.tag = element_text(size = 12)) ->

# --- Simulate Raw Stock Barcodes With Clustering -------------------------

## Saving for later
## stock_raw_barcode_ranking %>% filter(ranking <= 5000) %>% mutate(pct = count / sum(count))
## > stock_raw_barcode_ranking %>% filter(ranking <= 5000) %>% pull(count) %>% sum()

expand_L_dist <- function(kmer, L) {

    if (L == 0) {
        return(kmer)
        }

    original_kmer = kmer
    edit_dist_L = FALSE

    # while loop is to verify that edited kmer is actually L away from original kmer
    while (kmer == original_kmer) {
        for (i in 1:L) {
            nucs = c('A', 'C', 'G', 'T')
            length_kmer = str_length(kmer)
            ## print('kmer:', kmer)
            ## print(paste('length kmer:', length_kmer))

            edit = sample(c('insert', 'delete', 'sub'), 1)

            if (edit == 'insert') {
                ## print('case_insert')
                pos_insert = sample.int(n = length_kmer + 1, size = 1) - 1
                insert_nuc = sample(nucs, 1)
                kmer = paste0(str_sub(kmer, 0, pos_insert), insert_nuc, str_sub(kmer, pos_insert + 1))
                ## print(kmer)
            } else if (edit == 'delete') {
                ## print('case_delete')
                pos_del = sample.int(n = length_kmer, size = 1)
                kmer = paste0(str_sub(kmer, 0, pos_del-1), str_sub(kmer, pos_del+1))
                ## print(kmer)
            } else if (edit == 'sub') {
                ## print('case_sub')
                pos_sub = sample.int(n = length_kmer, size = 1)
                substring(kmer, pos_sub, pos_sub) = sample(nucs[! nucs %in% substring(kmer, pos_sub, pos_sub)], 1)
                # ^ Complicated expression inside sample is to avoid substituting a letter with itself
                ## print(kmer)
            }
            ## print(kmer)
        }
        ## edit_dist_L = ! (adist(original_kmer, kmer)[1, 1] == L)
    }
    return(kmer)
}

add_N_error_bcs_for_L <- function(barcode, N, L) {
    # The <<- is for assignment to global variables, which R discourages
    # But this seems like the most memory-efficient way of dealing with something
    # that probably requires a lot of memory

    ## if (N == 0) {return(NULL)}

    for (i in 1:N) {
        if (N == 0) {break}
        ## print(N)
        ## print(L)
        error_bc = expand_L_dist(barcode, L)
        
        if (error_bc %in% names(random_w_error_counts)) {
            ## print('in name exists')
            random_w_error_counts[error_bc] <<- random_w_error_counts[error_bc] + 1
            }
        else {
                ## print('name doesnot exist')
            random_w_error_counts[error_bc] <<- 1
            }
        }
    ## return(count_list)
    }

random_w_error_counts = c()

random_bc_counts = tibble(barcode = random_bcs,
                          count = top_n(stock_raw_counts, 5000, wt = count) %>%
                              arrange(desc(count)) %>%
                              pull(count)
                          )

random_bcs_w_errors_wo0 = one_plasmid_L_dist_pcts %>%
    select(distance, pct_of_1P) %>%
    filter(distance > 0) %>%
    expand_grid(., random_bc_counts) %>%
    select(barcode, count, distance, pct_of_1P) %>%  # these two lines purely for aesthetics
    arrange(desc(count)) %>%
    mutate(count_to_generate = floor(count * pct_of_1P))
    
## for (i in 1:nrow(head(random_bcs_w_errors_wo0)) {
##     print(random 
##     }

# Add error bcs to random_w_error_counts vector, exclude the original barcodes, as we can just add their numbers later.
for (i in 1:nrow(random_bcs_w_errors_wo0)) {
    add_N_error_bcs_for_L(random_bcs_w_errors_wo0$barcode[i],
                          random_bcs_w_errors_wo0$count_to_generate[i],
                          random_bcs_w_errors_wo0$distance[i])
    }

for (i in 1:nrow(random_bc_counts)) {
    if (random_bc_counts$barcode[[i]] %in% names(random_w_error_counts)) {
        random_w_error_counts[random_bc_counts$barcode[[i]]] = random_w_error_counts[random_bc_counts$barcode[[i]]] + random_bc_counts$count[[i]]
        }
    else {
        random_w_error_counts[random_bc_counts$barcode[[i]]] = random_bc_counts$count[[i]]
        }
    }

enframe(random_w_error_counts, name = 'barcode', value = 'count') %>%
    write_tsv('simulations/raw_Sep12.tsv', col_names = FALSE)

# Starcode was run, reading in results
L3_random_w_error_counts = read_tsv('simulations/L3_clustered_Sep12.tsv',
                                    col_names=c('barcode', 'count', 'elements')) %>%
    select(- elements)

L2_random_w_error_counts = read_tsv('simulations/L2_clustered_Sep12.tsv',
                                    col_names=c('barcode', 'count', 'elements')) %>%
    select(- elements)

# Count overlap of random 5k errors and L3 clustering of the generated errors
intersect(pull(L3_random_w_error_counts, barcode), pull(random_bc_counts, barcode)) %>% length()

intersect(pull(L2_random_w_error_counts, barcode), pull(random_bc_counts, barcode)) %>% length()

L3_cum99pct_random_w_error_counts = L3_random_w_error_counts %>%
    mutate(total = sum(count),
           pct = count / total,
           cum_pct = cumsum(pct)) %>%
    filter(cum_pct <= 0.99) 

intersect(pull(L3_cum99pct_random_w_error_counts, barcode), pull(random_bc_counts, barcode)) %>% length()

# Get one_plasmid_counts from the controls code, transferred via RDS
one_plasmid_counts = readRDS('one_plasmid_counts.rds')

ten_control_barcodes = c(BC7_1='TCACAGGGGTAA', BC7_2='ACAAGACCGGAA',
                         BC7_3='ATATAGAGCTGT', BC7_4='ACATACCTGCTA',
                         BC7_5='GTGTCAGGCACA', BC7_6='TGCCACTCTAGC',
                         BC7_8='CTCGATTCACTC', BC7_9='GAACCCGTGGAA',
                         BC7_12='CTGTATATTTTA', BC7_15='GAAACCATGACA')

raw_BC1_103_bcs = one_plasmid_counts %>%
    filter(cluster == 'Raw' & library == 'BC1_103') %>%
    pull(barcode)

L_dist_raw_BC1_103_BC7_1 = adist(ten_control_barcodes[['BC7_1']], raw_BC1_103_bcs)

one_plasmid_L_dist_pcts = one_plasmid_counts %>%
    filter(cluster == 'Raw' & library == 'BC1_103') %>%
    add_column(distance = as.vector(L_dist_raw_BC1_103_BC7_1)) %>%
    filter(distance <= 10) %>%
    group_by(distance) %>%
    summarize(count = sum(count)) %>%
    mutate(pct_of_1P = count / (filter(., distance == 0) %>% pull(count)))

# --- Experiment With Taking Top N Barcodes -------------------------------

plasmid_cluster_counts %>%
    filter(cluster == 'Raw' & barcode != '') %>%
    top_n(5000, count) %>%
    ## mutate(total = sum(count),
    ##        pct = count /total,
    ##        cum_pct = cumsum(pct)) %>%
    ## filter(cum_pct <= 0.99) %>%
    pull(barcode) %>%
    set_names(., nm=.) %>%
    adist() %>%
    as.data.frame() %>%
    rownames_to_column('stock_1') %>%
    as_tibble() %>%
    pivot_longer(-stock_1, names_to='stock_2', values_to='dist') %>%
    filter(stock_1 < stock_2) %>%
    ggplot(., aes(dist)) +
    geom_histogram(binwidth = 1, fill = 'goldenrod') +
    labs(x = 'Levenshtein Distance',
         y = 'Count of Pairwise Distances') +
    scale_x_continuous(breaks = seq(0, 12, 4)) +
    scale_y_continuous(breaks = NULL) +
        theme_minimal() -> p_top5k_plasmid_Raw_distance_distribution

ggsave('top5k_plasmid_Raw_distance_distribution.pdf',
       path='../plots/stock',
       plot = p_top5k_plasmid_Raw_distance_distribution)

## Virtually identical plot to L3 cumulative 99% plot


# --- Similarity Between Plasmid And Stock --------------------------------

# Shamelessly stolen from https://stats.stackexchange.com/questions/31565/compute-a-cosine-dissimilarity-matrix-in-r
# But modified to handle the incoming table has samples as columns and observations as rows
cosine_similarity <- function(table) {
    matrix = as.matrix(table)
    sim = t(t(matrix) / sqrt(colSums(matrix * matrix)))
    return(t(sim) %*% sim)
    }

plasmid_cluster_counts %>%
    filter(cluster == 'Raw' & barcode != '') %>%
    mutate(library = 'Plasmid') %>%
    select(library, barcode, count) %>%
    bind_rows(mutate(stock_raw_counts, library = 'Stock')) %>% 
    pivot_wider(names_from = library, values_from = count, values_fill = 0) %>%
    column_to_rownames(var = 'barcode') %>%
    cosine_similarity()

