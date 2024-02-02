library(patchwork)
library(ggridges)
library(ComplexUpset)
library(pheatmap)
library(corrplot)
library(tidyverse)

library(entropy)

library(kableExtra)

animal_colors = c('FL' = 'palevioletred3', 'FR' = 'palevioletred1', 'ML' = 'skyblue3', 'MR' = 'skyblue1')

# --- Get Urine Raw Counts ------------------------------------------------

JA19161_urine_raw_counts = list.files('../barcodes_e0.1/JA19161',
                                     pattern='D.*_cutadapt_counts.txt',
                                     full.names=TRUE) %>%
    set_names(., gsub('.*\\/(.*?)_.*', '\\1', ., perl=TRUE)) %>%
    imap(~ read_table(.x, col_names=c(paste0('JA19161_', .y), 'barcode'), col_types='ic')) %>%
    reduce(full_join, by='barcode')

JA19375_urine_raw_counts = list.files('../barcodes_e0.1/JA19375',
                                     pattern='D.*_cutadapt_counts.txt',
                                     full.names=TRUE) %>%
    set_names(., gsub('.*\\/(.*?)_.*', '\\1', ., perl=TRUE)) %>%
    imap(~ read_table(.x, col_names=c(paste0('JA19375_', .y), 'barcode'), col_types='ic')) %>%
    reduce(full_join, by='barcode')

urine_raw_counts = full_join(JA19161_urine_raw_counts, JA19375_urine_raw_counts) %>%
    pivot_longer(-barcode, names_to='run_sample', values_to='count') %>%
    replace_na(list(count = 0)) %>%      # replace NAs in count with 0s
    mutate(sample = gsub('JA.*_(.*?)(2|bis)?$', '\\1', run_sample)) %>%     # identify sample names like LUFR and LUFR2
    group_by(sample, barcode) %>%
    summarize(count = sum(count)) %>%
    ungroup()


# --- Get Tissue Raw Counts -----------------------------------------------

## Don't use PLMR! Some issue with that sample.
## For the pattern we want to exclude names that begin with D (urine), and also BC (stock)
JA19161_tissue_raw_counts = list.files('../barcodes_e0.1/JA19161',
                                      pattern='^[B-CE-Z][^C].*_cutadapt_counts.txt',
                                      full.names=TRUE) %>%
    set_names(., gsub('.*\\/(.*?)_.*', '\\1', ., perl=TRUE)) %>%
    imap(~ read_table(.x, col_names=c(paste0('JA19161_', .y), 'barcode'), col_types='ic')) %>%
    reduce(full_join, by='barcode') %>%
    select(-JA19161_PLMR)  # Bad particular sample

JA19375_tissue_raw_counts = list.files('../barcodes_e0.1/JA19375',
                                      pattern='^[B-CE-Z][^C].*_cutadapt_counts.txt',
                                      full.names=TRUE) %>%
    set_names(., gsub('.*\\/(.*?)_.*', '\\1', ., perl=TRUE)) %>%
    imap(~ read_table(.x, col_names=c(paste0('JA19375_', .y), 'barcode'), col_types='ic')) %>%
    reduce(full_join, by='barcode')

tissue_raw_counts = full_join(JA19161_tissue_raw_counts, JA19375_tissue_raw_counts) %>%
    pivot_longer(-barcode, names_to='run_sample', values_to='count') %>%
    replace_na(list(count = 0)) %>%      # replace NAs in count with 0s
    mutate(sample = gsub('JA.*_(.*?)(2|bis)?$', '\\1', run_sample)) %>%     # identify sample names like LUFR and LUFR2
    group_by(sample, barcode) %>%
    summarize(count = sum(count)) %>%
    ungroup()


# --- Read L3 Clustered Stock, Apply Cutoff -------------------------------

clustered_combined_stock_barcodes = read_tsv('../stock/BC7_combined_stock_mp_L3_clusters.tsv',
                                             col_names = c('barcode', 'count', 'elements'),
                                             progress = FALSE) %>%
    select(-elements)

cutoff_99pct_stock_barcodes = clustered_combined_stock_barcodes %>%
    mutate(total = sum(count),
           pct = count / total,
           cum_pct = cumsum(pct)) %>%
    filter(cum_pct <= 0.99)


# --- Compute Distances Between Sample And Stock Barcodes -----------------
## Cut off distances past threshold, choose smallest distances to stock barcode

cutoff_99pct_stock_barcodes %>%
    pull(barcode) %>%
    saveRDS('split_sample_barcodes/stock_cum99pct_cutoff.rds')

## Split unique urine barcodes into smaller pieces to parallelize adist computation
## This splitting factor splits into list of 10,000 barcode vectors

unique_sample_barcodes = c(urine_raw_counts$barcode, tissue_raw_counts$barcode) %>%
    unique()

sample_splitting_factor = ceiling(seq(1, length(unique_sample_barcodes)) / 10000)

split(unique_sample_barcodes, sample_splitting_factor) %>%
    imap(~ saveRDS(.x, paste0('split_sample_barcodes/sample_barcodes_', .y, '.rds')))

sym_diff <- function(a,b) setdiff(union(a,b), intersect(a,b))
## Used this to spot check a couple lines in the distance RDS files, to see if the stock barcodes are the same

# Distance to stock barcode threshold
threshold = 3

## Read in the table of sample barcodes to stock barcodes, with distances, restrict to min distance, then count
## number of stock barcodes equidistant to single sample barcode

sample_barcodes_within_threshold_distance = list.files('split_sample_barcodes',
                                                       pattern = '.*dist_[[:digit:]]+\\.rds',
                                                       full.names = TRUE) %>%
    map(~ readRDS(.) %>% filter(., dist <= threshold)) %>%
    map(~ group_by(., sample) %>% filter(., dist == min(dist))) %>%
    reduce(bind_rows) %>%
    mutate(n = n()) %>%
    rename(sample_bc = sample, stock_bc = stock)

## Associate stock barcode to every individual sample barcode, then divide sample count by the number of
## stock barcodes equidistance from the sample barcode. Suppose a sample barcode has count 3, and there are
## three stock barcodes equidistance (with minimum distance) from sample barcode; then 1 is assigned to each
## stock barcode from the original sample count of 3.

urine_to_stock_barcodes = inner_join(urine_raw_counts,
                                     sample_barcodes_within_threshold_distance,
                                     by = c('barcode' = 'sample_bc')) %>%
    mutate(count_div_n = count / n) %>%
    select(sample, stock_bc, count_div_n) %>%
    group_by(sample, stock_bc) %>%
    summarize(weighted_count = sum(count_div_n)) # Sum the weighted values of the same stock barcodes

    ## ^ In above code block, many distinct sample barcodes were converted to the same stock barcodes,
    ## so we need to sum those up, which the group_by and mutate do.

tissue_to_stock_barcodes = inner_join(tissue_raw_counts,
                                     sample_barcodes_within_threshold_distance,
                                     by = c('barcode' = 'sample_bc')) %>%
    mutate(count_div_n = count / n) %>%
    select(sample, stock_bc, count_div_n) %>%
    group_by(sample, stock_bc) %>%
    summarize(weighted_count = sum(count_div_n)) # Sum the weighted values of the same stock barcodes


# --- Urine Analysis ------------------------------------------------------
    
## Read in data on concentrations of barcodes in urine,
## calculate barcode levels based on raw counts and concentration per sample

urine_pcr = read_tsv('26-04-2020_qPCR_data_urines.tsv', skip=7) %>%
    select(-starts_with('X')) %>%
    select(-starts_with('Copies')) %>%
    select(-Comments) %>%
    rename(Days_pi = `DAYS p.i.`,
           animal = Animals,
           sample = `Sample Name`,
           mean_ul_urine = `Quantity Mean/µl urine`,
           genomes_reaction = `Genomes eq copies/ enrich PCR reaction`,
           genomes_ul_urine = `Genomes eq copies/ µl urine`,
           quantity_mean = `Quantity Mean`,
           DNA_vs_urine_conc = `DNA conc vs urine conc factor`,
           vol_DNA_pcr = `Vol DNA (µl) used in PCR`,
           vol_urine_pcr = `Vol eq urine (µl) used in PCR`
           ) %>%
    select(sample, animal, Days_pi, mean_ul_urine, everything())
           
urine_bc_levels = urine_to_stock_barcodes %>%
    group_by(sample) %>%
    mutate(frac_weighted_count = weighted_count / sum(weighted_count)) %>%
    ungroup() %>%
    left_join(urine_pcr) %>%
    mutate(bc_level = frac_weighted_count * mean_ul_urine) %>%
    select(sample, animal, Days_pi, stock_bc, bc_level, frac_weighted_count, mean_ul_urine) %>%
    rename(barcode = stock_bc)

## -- Common Functions ----------------------------------------------------

top_barcodes_by_max <- function(barcode_table, top_n) {
    top_barcode_levels = barcode_table %>%
        group_by(animal, barcode) %>%      # Get highest level for each barcode across days
        mutate(max_level = max(bc_level)) %>%
        ungroup() %>%
        ## select(-frac_weighted_count) %>%
        nest_by(animal, max_level) %>%     # Nest all the longitudinal data to focus on max levels
        group_by(animal) %>%               # Get rid of rowwise grouping
        slice_max(max_level, n = top_n) %>%   # Get max 10 bc levels per animal
        ungroup() %>% 
        unnest(cols = data)                # Expand all the data that corresponds to just the max levels
    return(top_barcode_levels)
    }

top_barcodes_by_frac_max <- function(barcode_table, top_n) {
    top_barcode_levels = barcode_table %>%
        group_by(animal, barcode) %>%      # Get highest level for each barcode across days
        mutate(max_level = max(frac_weighted_count)) %>%
        ungroup() %>%
        ## select(-frac_weighted_count) %>%
        nest_by(animal, max_level) %>%     # Nest all the longitudinal data to focus on max levels
        group_by(animal) %>%               # Get rid of rowwise grouping
        slice_max(max_level, n = top_n) %>%   # Get max 10 bc levels per animal
        ungroup() %>% 
        unnest(cols = data)                # Expand all the data that corresponds to just the max levels
    return(top_barcode_levels)
    }

slice_barcodes_by_max <- function(barcode_table, top_rank, bottom_rank) {
    top_barcode_levels = barcode_table %>%
        group_by(animal, barcode) %>%      # Get highest level for each barcode across days
        mutate(max_level = max(bc_level)) %>%
        ungroup() %>%
        ## select(-frac_weighted_count) %>%
        nest_by(animal, barcode, max_level) %>%     # Nest all the longitudinal data to focus on max levels
        group_by(animal) %>%               # Get rid of rowwise grouping
        arrange(desc(max_level)) %>%
        slice(top_rank:bottom_rank) %>%   # Get max 10 bc levels per animal
        ungroup() %>% 
        unnest(cols = data)                # Expand all the data that corresponds to just the max levels
    return(top_barcode_levels)
    }

## -- Urine Area Plots Of Top 10 Winners ----------------------------------

individual_winner_area_plots = urine_bc_levels %>%
    top_barcodes_by_max(10) %>%
    split(.$animal) %>%
    imap(~ mutate(., barcode = fct_reorder(as.factor(barcode), max_level, .desc=TRUE))) %>%
    imap(~ ggplot(.x, aes(x = Days_pi, fill = barcode)) +
             geom_area(aes(y = mean_ul_urine), fill = 'grey85') +
             geom_area(aes(y = bc_level)) +
             facet_wrap(vars(barcode), ncol = 1) +
             labs(title = 'Top 10 Winners',
                  subtitle = paste('Mouse', .y),
                  x = 'Days post-injection',
                  y = 'Barcode level') +
             theme_minimal() +
             theme(legend.position = 'none'))

summed_winner_area_plots = urine_bc_levels %>%
    top_barcodes_by_max(10) %>%
    group_by(animal, Days_pi, mean_ul_urine) %>%
    summarize(sum_winner_barcodes = sum(bc_level)) %>%
    split(.$animal) %>%
    imap(~ ggplot(.x, aes(x = Days_pi)) +
             geom_area(aes(y = mean_ul_urine), fill = 'grey85') +
             geom_area(aes(y = sum_winner_barcodes), fill = 'navajowhite3') +
             labs(title = 'Sum Of Top 10 Winner Barcodes',
                  x = 'Days post-injection',
                  y = 'Barcode level') +
             theme_minimal() +
             theme(legend.position = 'none'))

pwalk(list(names(individual_winner_area_plots), individual_winner_area_plots, summed_winner_area_plots),
      ~ ggsave(paste0('urine_winners_w_sum_linear_', ..1, '_area.pdf'),
               plot = ..2 / ..3 +
                   plot_layout(heights = c(10, 1)) +
                   plot_annotation(
                       title = paste('Top 10 Barcodes In Mouse', ..1),
                       tag_levels = 'A',
                       theme = theme(plot.title = element_text(size=18))
                   ),
               path = '../plots/urine_tissue_manuscript',
               height = 14)
      )

## -- Number Of Unique Barcodes -------------------------------------------

num_unique_barcode_plots = urine_bc_levels %>%
    filter(bc_level > 0) %>%
    group_by(animal, Days_pi) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    split(.$animal) %>%
    imap(~ ggplot(.x, aes(x = Days_pi, y = n)) +
             geom_line(color = animal_colors[[.y]]) +
             labs(title = 'Number Of Unique Barcodes',
                  x = 'Days post-injection',
                  y = 'Count') +
             theme_minimal()
         ) %>%
    iwalk(~ ggsave(paste0('urine_num_unique_bcs_', .y, '_line.pdf'),
                  plot = .x,
                  path = '../plots/urine_tissue_manuscript',
                  height = 4)
         )

urine_bc_levels %>%
    filter(bc_level > 0) %>%
    group_by(animal, Days_pi) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    ggplot(aes(x = Days_pi, y = n)) +
    geom_line(aes(color = as_factor(animal)),
              linewidth = 0.8,
              show.legend = FALSE) +
    scale_color_manual(values = animal_colors) +
    facet_wrap(vars(animal), ncol = 1) +
             labs(title = 'Number Of Unique Barcodes',
                  x = 'Days post-injection',
                  y = 'Count') +
             theme_minimal()

ggsave('urine_num_unique_bcs_all_animals.pdf',
       path = '../plots/urine_tissue_manuscript')

    iwalk(~ ggsave(paste0('urine_num_unique_bcs_', .y, '_line.pdf'),
                  plot = .x,
                  path = '../plots/urine_tissue_manuscript',
                  height = 4)
         )
    
    
pwalk(list(names(individual_winner_area_plots), individual_winner_area_plots, summed_winner_area_plots, num_unique_barcode_plots),
      ~ ggsave(paste0('urine_winners_w_sum_linear_and_n_unique_', ..1, '_area.pdf'),
               plot = ..2 / ..3 / ..4 +
                   plot_layout(heights = c(10, 1, 1)) +
                   plot_annotation(
                       title = paste('Top 10 Barcodes In Mouse', ..1),
                       tag_levels = 'A',
                       theme = theme(plot.title = element_text(size=18))
                   ),
               path = '../plots/urine_and_tissue',
               height = 16)
      )

urine_bc_levels %>%
    split(.$animal) %>%
    imap(~ group_by(.x, barcode)) %>%
    imap(~ summarize(.x, mean = mean(bc_level))) %>%
    imap(~ ggplot(.x, aes(sqrt(mean))) +
             geom_histogram(bins = 100)) %>%
    iwalk(~ ggsave(paste0('mean_bc_level_per_barcode_', .y, '_hist.pdf'),
                   path = '../plots/urine_and_tissue'))

# Numbered time points
urine_bc_levels %>%
    filter(animal == 'FL') %>%
    distinct(Days_pi) %>%
    arrange(Days_pi) %>%
    pull(Days_pi) %>%
    set_names(1:length(.))

num_timepoints_per_animal = urine_bc_levels %>%
    filter(bc_level > 0) %>%
    select(animal, Days_pi) %>%
    distinct() %>%
    group_by(animal) %>%
    summarize(n = n()) %>%
    deframe()

# Getting tibble of list of unique barcodes
nested_barcodes = urine_bc_levels %>%
    filter(animal == 'FL' & bc_level > 0) %>%
    select(animal, Days_pi, barcode) %>%
    nest_by(animal, Days_pi) %>%
    mutate(data = map(data, as_vector)) 

timepoint_intersections = expand_grid(first = 1:num_timepoints_per_animal[['FL']],
                                      second = 1:num_timepoints_per_animal[['FL']]) %>%
    filter(first < second) %>%
    mutate(distance = second - first) %>%
    rowwise() %>%
    mutate(unique_bcs_in_intersection = length(intersect(
               nested_barcodes[[first, 'data']][[1]],
               nested_barcodes[[second, 'data']][[1]]
               ))
           )

timepoint_intersections %>%
    ggplot(aes(x = factor(distance), y = unique_bcs_in_intersection)) +
    geom_boxplot() +
    theme_minimal() 

    ggsave('test_intersections_FL_box.pdf',
           path = '../plots/urine_and_tissue',
           width = 10)
                      
    ggsave('test_intersections_FL_violin.pdf',
           path = '../plots/urine_and_tissue')
                      
    ggsave('test_intersections_FL_dot.pdf',
           path = '../plots/urine_and_tissue')

timepoint_interesections = num_timepoints_per_animal %>%
    imap(~ expand_grid(first = 1:.x, second = 1:.x)) %>% map(~ filter(.x, first < second))

loser_400_nested_barcodes = urine_bc_levels %>%
    filter(animal == 'FL' & bc_level > 0) %>%
    group_by(Days_pi) %>%
    slice_min(order_by = bc_level, n = 400) %>%
    ungroup() %>%
    select(animal, Days_pi, barcode) %>%
    nest_by(animal, Days_pi) %>%
    mutate(data = map(data, as_vector))
    
loser_400_timepoint_intersections = expand_grid(first = 1:num_timepoints_per_animal[['FL']],
                                                second = 1:num_timepoints_per_animal[['FL']]) %>%
    filter(first < second) %>%
    mutate(distance = second - first) %>%
    rowwise() %>%
    mutate(unique_bcs_in_intersection = length(intersect(
               loser_400_nested_barcodes[[first, 'data']][[1]],
               loser_400_nested_barcodes[[second, 'data']][[1]]
               ))
           )

loser_400_timepoint_intersections %>%
    ggplot(aes(x = distance, y = unique_bcs_in_intersection)) +
    geom_point() +
    geom_smooth() +
    ## ylim(0, 400) +
    theme_minimal() 

ggsave('loser_400_intersections_FL_box.pdf',
        path = '../plots/urine_and_tissue',
        width = 10)

ggsave('loser_400_intersections_nolimits_FL_dot_loess.pdf',
        path = '../plots/urine_and_tissue',
        width = 10)

animal_days = urine_bc_levels %>%
    select(animal, Days_pi) %>%
    distinct() %>%
    split(.$animal) %>%
    imap(~ arrange(.x, Days_pi)) %>%
    imap(~ pull(.x, Days_pi))

loser_400_intersections_day_distance = expand_grid(first = animal_days[['FL']],
                                                   second = animal_days[['FL']]) %>%
    filter(first < second) %>%
    mutate(distance = second - first) %>%
    rowwise() %>%
    mutate(unique_bcs_in_intersection = length(
               intersect(filter(loser_400_nested_barcodes, Days_pi == first)$data[[1]],
                         filter(loser_400_nested_barcodes, Days_pi == second)$data[[1]]
                         )
           )
           )

loser_400_intersections_day_distance %>%
    ggplot(aes(x = distance, y = unique_bcs_in_intersection)) +
    geom_point() +
    geom_smooth() +
    ylim(0, 400) +
    theme_minimal() 

ggsave('loser_400_intersections_days_FL_dot_loess.pdf',
        path = '../plots/urine_and_tissue',
        width = 10)

# Now for all the animals
loser_400_nested_barcodes = urine_bc_levels %>%
    filter(bc_level > 0) %>%
    group_by(animal, Days_pi) %>%
    slice_min(order_by = bc_level, n = 400) %>%
    ungroup() %>%
    select(animal, Days_pi, barcode) %>%
    nest_by(animal, Days_pi) %>%
    mutate(data = map(data, as_vector))
    
loser_400_intersections_day_distance = animal_days %>%
    imap(~ expand_grid(first = animal_days[[.y]],
                       second = animal_days[[.y]])
         ) %>%
    imap(~ filter(.x, first < second)) %>%
    imap(~ mutate(.x, distance = second - first)) %>%
    imap(~ rowwise(.x)) %>%
    imap(~ mutate(.x, unique_bcs_in_intersection = length(
                          intersect(
                              filter(loser_400_nested_barcodes,
                                     animal == .y & Days_pi == first)$data[[1]],
                              filter(loser_400_nested_barcodes,
                                     animal == .y & Days_pi == second)$data[[1]]
                          )
                      )
                  )
         )

loser_400_intersections_day_distance %>%
    imap(~ ggplot(.x, aes(x = distance, y = unique_bcs_in_intersection)) +
             geom_point(color = animal_colors[[.y]]) +
             geom_smooth(color = 'goldenrod') +
             ## ylim(0, 400) +
             theme_minimal()
         ) %>%
    imap(~ ggsave(paste('loser_400_intersections_days_nolimits', .y, 'dot_loess.pdf', sep = '_'),
                  plot = .x,
                  path = '../plots/urine_and_tissue',
                  width = 10)
         )

# Experiment With Plot For Checking Barcodes Bouncing Above And Below Median

FL_barcodes_below_median = urine_bc_levels %>%
    filter(animal == 'FL' & bc_level <= 8.174662e-05) %>%
    select(animal, Days_pi, barcode) %>%
    group_by(animal, Days_pi) %>%
    nest() %>%
    mutate(data = map(data, as_vector)) %>%
    ungroup() %>%
    select(-animal)

%>%
    deframe()
    

expand_grid(first = animal_days[['FL']], second = animal_days[['FL']]) %>%
    filter(first < second) %>%
    rowwise() %>%
    mutate(first_bcs = filter(FL_barcodes_below_median, Days_pi == first)$data) %>%
    mutate(second_bcs = filter(FL_barcodes_below_median, Days_pi == second)$data)

%>%
    mutate(intersection = intersect(first_bcs, second_bcs))

FL_successive_below_median = tibble(
    day = animal_days[['FL']][2:length(animal_days[['FL']])],
    prev_day = animal_days[['FL']][1:length(animal_days[['FL']])-1]
) %>%
    left_join(FL_barcodes_below_median, by = c('day' = 'Days_pi')) %>%
    left_join(FL_barcodes_below_median, by = c('prev_day' = 'Days_pi')) %>%
    mutate(successive_intersection = FL_successive_intersections[[2:length(FL_successive_intersections)]])
              
    rowwise() %>%
    mutate(num_below_median = as.character(FL_barcodes_below_median[[day]]))
    
    nest_join(FL_barcodes_below_median, by = c('second' = 'Days_pi'), name = 'second_bcs') %>% pull(second_bcs) %>% pluck(1) 

# Count successive intersection of bcs below median
FL_successive_intersections = FL_barcodes_below_median %>% arrange(Days_pi) %>% pull(data) %>% accumulate(intersect) %>% map(~ length(.)) %>% flatten_int()

FL_successive_intersections %>%
    enframe() %>%
    mutate(Days_pi = urine_bc_levels %>% filter(animal == 'FL') %>% distinct(Days_pi) %>% arrange(Days_pi) %>% pull()) %>%
    select(-name) %>%
    ggplot(aes(x = Days_pi, y = value)) +
    geom_line() +
    labs(title = 'Successive Intersections Of Barcodes Below Median Level In Mouse FL',
         x = 'Days post-injection',
         y = 'Number of barcodes') +
    theme_minimal()

ggsave('successive_intersections_FL.pdf',
       path = '../plots/urine_and_tissue')

# Intersections of barcodes above and below threshold

threshold = urine_bc_levels %>% filter(animal == 'FL') %>% pull(bc_level) %>% median()

FL_barcodes_below_threshold = urine_bc_levels %>%
    filter(animal == 'FL' & bc_level > 0 & bc_level <= threshold) %>%
    select(animal, Days_pi, barcode)

%>%
    group_by(animal, Days_pi) %>%
    nest() %>%
    mutate(data = map(data, as_vector)) %>%
    ungroup() %>%
    select(-animal) %>%
    arrange(Days_pi) %>%
    deframe()

FL_barcodes_above_threshold = urine_bc_levels %>%
    filter(animal == 'FL' & bc_level > threshold) %>%
    select(animal, Days_pi, barcode) %>%
    group_by(animal, Days_pi) %>%
    nest() %>%
    mutate(data = map(data, as_vector)) %>%
    ungroup() %>%
    select(-animal) %>%
    arrange(Days_pi) %>%
    deframe()

intersections_of_below_thresholds = map2(FL_barcodes_below_threshold[2:length(FL_barcodes_below_threshold)],
                                         FL_barcodes_below_threshold[1:length(FL_barcodes_below_threshold)-1],
                                         ~intersect(.x, .y)) %>%
    map(~ length(.))
    
intersections_of_above_below_thresholds = map2(FL_barcodes_below_threshold[2:length(FL_barcodes_below_threshold)],
                                               FL_barcodes_above_threshold[1:length(FL_barcodes_above_threshold)-1],
                                               ~intersect(.x, .y)) %>%
    map(~ length(.))
    
## intersections_of_above_below_thresholds = 

# Check min barcode levels relative to median

urine_bc_levels %>%
    filter(animal == 'FL') %>%
    filter(bc_level > 0) %>%
    group_by(Days_pi) %>%
    summarize(min_level = min(bc_level)) %>%
    ggplot(aes(x = Days_pi, y = min_level)) +
    geom_line(color = animal_colors[['FL']]) +
    geom_hline(yintercept = threshold, color = 'gray60') +
    theme_minimal()

ggsave('min_bc_level_FL_line.pdf',
       path = '../plots/urine_and_tissue')


## -- Compare Stock And Day 1 Repertoires ---------------------------------

# Shamelessly stolen from https://stats.stackexchange.com/questions/31565/compute-a-cosine-dissimilarity-matrix-in-r
# But modified to handle the incoming table has samples as columns and observations as rows
cosine_similarity <- function(table) {
    matrix = as.matrix(table)
    sim = t(t(matrix) / sqrt(colSums(matrix * matrix)))
    return(t(sim) %*% sim)
    }

first_day_and_stock_frac_counts = urine_bc_levels %>%
    filter(Days_pi == 1) %>%
    select(animal, barcode, frac_weighted_count) %>%
    rename(repertoire = animal, pct = frac_weighted_count) %>%
    bind_rows(
        cutoff_99pct_stock_barcodes %>%
        mutate(total = sum(count), pct = count / total) %>%
        select(barcode, pct) %>%
        mutate(repertoire = 'stock')
    ) %>%
    pivot_wider(names_from = repertoire, values_from = pct, values_fill = 0) %>%
    column_to_rownames(var = 'barcode')

pdf('../plots/urine_tissue_manuscript/day1_and_stock_cosine_corrplot.pdf', onefile = FALSE)
    cosine_similarity(first_day_and_stock_frac_counts) %>%
        corrplot.mixed(.,
                       upper = 'circle',
                       lower.col = 'grey',
                       tl.col = 'black',
                       is.corr = FALSE,
                       col.lim = c(0, max(.)),  # For some reason, explicitly setting upper limit to 1 gave error
                       upper.col = COL1('Purples', 200))
dev.off()

## -- Cosine Similarities Per Animal Across Time --------------------------

cos_sim_across_days = urine_bc_levels %>%
    filter(animal == 'FL') %>%
    select(Days_pi, barcode, frac_weighted_count) %>%
    rename(pct = frac_weighted_count) %>%
    arrange(Days_pi) %>%
    pivot_wider(names_from = Days_pi, values_from = pct, values_fill = 0) %>%
    column_to_rownames(var = 'barcode') %>%
    cosine_similarity(.)

days_past_first = sort(unique(urine_bc_levels$Days_pi))[-1]

scale_cos_sim_to_mean_ul_urine <- function(table, animal) {
    table = table %>% filter(animal == animal)
    scale_factor = ceiling(max(table$mean_ul_urine) / 1000) * 1000  # Multiples of 1000 seem appropriate for this particular data

    return(scale_factor)
    }
    
cos_sim_subsequent_timepoints = urine_bc_levels %>%
    ## filter(animal == 'FL') %>%
    select(animal, Days_pi, mean_ul_urine) %>%
    group_by(animal) %>%
    distinct() %>%
    arrange(Days_pi) %>%
    mutate(second_timepoint = lead(Days_pi)) %>%
    slice(1:n() - 1) %>%
    rowwise() %>%
    mutate(cos_sim = cos_sim_across_days[as.character(Days_pi), as.character(second_timepoint)])

cos_sim_across_days = urine_bc_levels %>%
    ## filter(animal == 'FL') %>%
    select(animal, Days_pi, barcode, frac_weighted_count) %>%
    rename(pct = frac_weighted_count) %>%
    split(.$animal) %>%
    imap(~ select(.x, -animal)) %>%
    imap(~ arrange(.x, Days_pi)) %>%
    imap(~ pivot_wider(.x, names_from = Days_pi, values_from = pct, values_fill = 0)) %>%
    imap(~ column_to_rownames(.x, var = 'barcode')) %>%
    imap(~ cosine_similarity(.x))

cos_sim_subsequent_timepoints = urine_bc_levels %>%
    select(animal, Days_pi, mean_ul_urine) %>%
    split(.$animal) %>%
    imap(~ distinct(.x)) %>%
    imap(~ arrange(.x, Days_pi)) %>%
    imap(~ mutate(.x, second_timepoint = lead(Days_pi))) %>%
    imap(~ slice(.x, 1:n() -1)) %>%
    imap(~ rowwise(.x)) %>%
    imap(~ mutate(.x, cos_sim = cos_sim_across_days[[.y]][as.character(Days_pi), as.character(second_timepoint)]))

cos_sim_subsequent_timepoints %>%
    imap(~ ggplot(.x, aes(x = second_timepoint)) +
            geom_area(aes(y = mean_ul_urine), stat = 'identity', fill = 'gray85') +
             geom_line(aes(y = cos_sim * scale_cos_sim_to_mean_ul_urine(urine_bc_levels, .y)),
                       color = animal_colors[.y]) +
            scale_x_continuous(name = 'Days Post-Injection') +
            scale_y_continuous(name = 'Genomes Per uL',
                            sec.axis = sec_axis(~ ./scale_cos_sim_to_mean_ul_urine(urine_bc_levels, .y),
                                                name = 'Cosine Similarity')
                            ) +
            labs(title = 'Cosine Similarity Between Subsequent Time Points vs Total Shed Virus',
                subtitle = paste('Mouse', .y)) +
             theme_minimal()
         ) %>%
    iwalk(~ ggsave(paste('cos_sim_subsequent_days', .y, 'line.pdf', sep = '_'),
                   path = '../plots/urine_and_tissue',
                   plot = .x,
                   height = 4,
                   width = 10)
          )

cos_sim_subsequent_timepoints %>%
    bind_rows() %>%
    ggplot(aes(x = second_timepoint)) +
    geom_area(aes(y = mean_ul_urine), stat = 'identity', fill = 'gray85') +
    geom_line(aes(y = cos_sim * scale_cos_sim_to_mean_ul_urine(urine_bc_levels),
                  color = as_factor(animal)),
              linewidth = 0.8,
              show.legend = FALSE) +
    scale_x_continuous(name = 'Days Post-Injection') +
    scale_y_continuous(name = 'Genomes Per uL',
                       sec.axis = sec_axis(~ ./scale_cos_sim_to_mean_ul_urine(urine_bc_levels, .y), name = 'Cosine Similarity')) +
    scale_color_manual(values = animal_colors) +
    facet_wrap(vars(animal), ncol = 1) +
    labs(title = 'Cosine Similarity Between Subsequent Time Points vs Total Shed Virus') +
    theme_minimal()

ggsave('cos_sim_subsequent_days_all_line.pdf',
       path = '../plots/urine_tissue_manuscript')

## -- Urine Ridge Plots ---------------------------------------------------

winner_barcode_levels = urine_bc_levels %>%
    group_by(animal, barcode) %>%      # Get highest level for each barcode across days
    mutate(max_level = max(bc_level)) %>%
    ungroup() %>%
    ## select(-frac_weighted_count) %>%
    nest_by(animal, max_level) %>%     # Nest all the longitudinal data to focus on max levels
    group_by(animal) %>%               # Get rid of rowwise grouping
    slice_max(max_level, n = 10) %>%   # Get max 10 bc levels per animal
    ungroup() %>% 
    unnest(cols = data)                # Expand all the data that corresponds to just the max levels

urine_bc_levels %>%
    top_barcodes_by_max(30) %>%
    split(.$animal) %>%
    imap(~ ggplot(.x, aes(x = Days_pi, y = barcode, height = bc_level)) +
             geom_ridgeline(fill = animal_colors[[.y]],
                            color = 'grey90',
                            scale = 0.0005) +
             scale_x_continuous(expand=c(0.01, 0)) +
             scale_y_discrete(expand=c(0.01, 0)) +
             labs(x = 'Days post-injection',
                  y = 'Barcode',
                  title = 'Top 30 Barcodes Ridge Plot',
                  subtitle = paste('Mouse', .y)) +
             theme_ridges()) %>%
    iwalk(~ ggsave(paste0('urine_top30_linear_', .y, '_ridge.pdf'),
                   plot = .x,
                   path = '../plots/urine_tissue_manuscript',
                   height = 17))
             
### - Using Ratio Of Counts Instead Of Barcode Level ----------------------

urine_bc_levels %>%
    top_barcodes_by_frac_max(30) %>%
    split(.$animal) %>%
    imap(~ ggplot(.x, aes(x = Days_pi, y = barcode, height = frac_weighted_count)) +
             geom_ridgeline(fill = animal_colors[[.y]],
                            color = 'grey90') +
             scale_x_continuous(expand=c(0.01, 0)) +
             scale_y_discrete(expand=c(0.01, 0)) +
             labs(x = 'Days post-injection',
                  y = 'Barcode',
                  title = 'Percent Of Total Counts, Top 30 Barcodes Selected By Percentage',
                  subtitle = paste('Mouse', .y)) +
             theme_ridges()) %>%
    iwalk(~ ggsave(paste0('urine_pct_top30_linear_', .y, '_ridge.pdf'),
                   plot = .x,
                   path = '../plots/urine_and_tissue',
                   height = 17))
             
             
## -- GC Plot -------------------------------------------------------------

tibble(animal = c('FL', 'FR', 'ML', 'MR'),        # Adding stock barcodes as timepoint 0
       Days_pi = c(0, 0, 0, 0)) %>%
    mutate(stock = map(.x = animal, ~ cutoff_99pct_stock_barcodes)) %>%
    unnest(., stock) %>%
    select(animal, Days_pi, barcode, pct) %>%
    bind_rows(urine_bc_levels %>%
              select(animal, Days_pi, barcode, frac_weighted_count) %>%
              rename(pct = frac_weighted_count)) %>%
    mutate(GC_pct = str_count(barcode, '[GC]') / str_length(barcode)) %>%
    mutate(wt_GC_pct = GC_pct * pct) %>%          # multiply %GC of barcode by %barcodes in sample
    group_by(animal, Days_pi) %>%
    summarize(timepoint_GC_pct = sum(wt_GC_pct)) %>%
    ggplot(aes(x = Days_pi, y = timepoint_GC_pct)) +
    geom_line(aes(color = as_factor(animal)),
              linewidth = 0.8,
              show.legend = FALSE) +
    facet_wrap(vars(animal), ncol = 1) +
    labs(title = 'Mean GC Content',
         x = 'Days Post-Injection',
         y = 'GC Content') +
    scale_color_manual(values = animal_colors) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format()) +
    theme_minimal()

ggsave('GC_content_line.pdf',
        path = '../plots/urine_tissue_manuscript')
    

### GC Content In Top 10 Barodes Per Mouse

urine_bc_levels %>%
    top_barcodes_by_max(10) %>%
    distinct(animal, barcode) %>%
    mutate(GC_pct = str_count(barcode, '[GC]') / str_length(barcode)) %>%
    select(animal, barcode, GC_pct) %>%
    ggplot(aes(animal, GC_pct)) +
    geom_boxplot(aes(color = animal), linewidth = 0.8, show.legend = FALSE) +
    labs(title = 'GC Content For Top 10 Barcodes In Urine',
         x = 'Mouse',
         y = 'GC Content') +
    scale_color_manual(values = animal_colors) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format()) +
    theme_minimal()

ggsave('GC_urine_winners_boxplot.pdf',
        path = '../plots/urine_tissue_manuscript')

gc_pct_all_barcodes = urine_bc_levels %>%
    distinct(barcode) %>%
    mutate(GC_pct = str_count(barcode, '[GC]') / str_length(barcode)) %>%
    add_column(animal = 'All Barcodes', .before = 1)

urine_bc_levels %>%
    top_barcodes_by_max(10) %>%
    distinct(animal, barcode) %>%
    mutate(GC_pct = str_count(barcode, '[GC]') / str_length(barcode)) %>%
    bind_rows(., gc_pct_all_barcodes) %>%
    mutate(animal = fct_relevel(animal, c('FL', 'FR', 'ML', 'MR', 'All Barcodes'))) %>%
    ggplot(aes(animal, GC_pct)) +
    geom_boxplot(aes(color = animal), linewidth = 0.8, show.legend = FALSE) +
    labs(title = 'GC Content For Top 10 Barcodes Per Animal In Urine',
         x = 'Mouse',
         y = 'GC Content') +
    scale_color_manual(values = c(animal_colors, `All Barcodes` = 'gray70')) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format()) +
    theme_minimal()

ggsave('GC_urine_winners_boxplot.pdf',
        path = '../plots/urine_tissue_manuscript')

tissue_bc_levels %>%
    top_barcodes_by_max(10) %>%
    distinct(animal, barcode) %>%
    mutate(GC_pct = str_count(barcode, '[GC]') / str_length(barcode)) %>%
    select(animal, barcode, GC_pct) %>%
    ggplot(aes(animal, GC_pct)) +
    geom_boxplot(aes(color = animal), linewidth = 0.8, show.legend = FALSE) +
    labs(title = 'GC Content For Top 10 Barcodes In Tissue',
         x = 'Mouse',
         y = 'GC Content') +
    scale_color_manual(values = animal_colors) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format()) +
    theme_minimal()

ggsave('GC_tissue_winners_boxplot.pdf',
        path = '../plots/urine_tissue_manuscript')


### Length Of Top 10 Barcodes

urine_bc_levels %>%
    top_barcodes_by_max(10) %>%
    distinct(animal, barcode) %>%
    mutate(length = str_length(barcode)) %>%
    select(animal, barcode, length) %>%
    ggplot(aes(animal, length)) +
    geom_jitter(aes(color = animal),
                height = 0,
                width = 0.25,
                show.legend = FALSE) +
    labs(title = 'Length Of Top 10 Barcodes In Urine',
         x = 'Mouse',
         y = 'Length') +
    scale_color_manual(values = animal_colors) +
    theme_minimal()

ggsave('length_urine_winners_jitter.pdf',
       path = '../plots/urine_tissue_manuscript',
       height = 5)

urine_bc_levels %>%
    top_barcodes_by_max(10) %>%
    distinct(animal, barcode) %>%
    mutate(length = str_length(barcode)) %>%
    select(animal, barcode, length) %>%
    ggplot(aes(length)) +
    geom_histogram(aes(fill = animal),
                   binwidth = 1,
                   breaks = seq(4, 12),
                   show.legend = FALSE) +
    facet_wrap(vars(animal), ncol = 1) +
    labs(title = 'Length Of Top 10 Barcodes In Urine',
         x = 'Mouse',
         y = 'Length') +
    scale_fill_manual(values = animal_colors) +
    theme_minimal()

ggsave('length_urine_winners_histogram.pdf',
       path = '../plots/urine_tissue_manuscript',
       height = 5)

urine_bc_levels %>%
    top_barcodes_by_max(10) %>%
    distinct(animal, barcode) %>%
    mutate(length = str_length(barcode)) %>%
    select(animal, barcode, length) %>%
    ## count(animal, length) %>%
    ggplot(aes(as_factor(length))) +
    geom_bar(aes(fill = animal),
                   ## binwidth = 1,
                   ## breaks = seq(4, 12),
                   show.legend = FALSE) +
    facet_wrap(vars(animal), ncol = 1) +
    labs(title = 'Length Of Top 10 Barcodes In Urine',
         x = 'Mouse',
         y = 'Length') +
    ## ylim(0, 10) +
    scale_y_discrete(limits = seq(0, 10, 2)) +
    scale_fill_manual(values = animal_colors) +
    theme_minimal()

ggsave('length_urine_winners_histogram.pdf',
       path = '../plots/urine_tissue_manuscript')



## -- Rank Urine Winners In Stock -----------------------------------------

urine_bc_levels %>%
    top_barcodes_by_max(10) %>%
    filter(bc_level == max_level) %>%
    split(.$animal) %>%
    imap(~ mutate(.x, urine_rank = desc(max_level) %>% min_rank() %>% as_factor())) %>%
    bind_rows() %>% 
    left_join(cutoff_99pct_stock_barcodes %>%
              mutate(stock_rank = desc(count) %>% min_rank())) %>%
    select(animal, max_level, barcode, urine_rank, count, stock_rank) %>%
    ggplot(aes(x = fct_relevel(urine_rank, rev), y = stock_rank)) +
    geom_point(aes(color = as_factor(animal)), size=3, show.legend=FALSE) +
    facet_wrap(vars(animal)) +
    labs(title = 'Rank Of Urine Winners In Stock',
         x = "Rank Of Winner In Urine",
         y = "Rank In Stock") +
    scale_y_continuous(limits = c(1, dim(cutoff_99pct_stock_barcodes)[[1]])) +
    scale_x_discrete(breaks = seq(1, 10)) +
    scale_color_manual(values=animal_colors) +
    coord_flip() +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(linetype = 'dotted',
                                            color = 'gray70'),
          panel.grid.minor.y = element_blank(),
          panel.spacing = unit(2, 'lines'),
          strip.background = element_rect(fill = 'gray90',
                                          color = 'gray90')
          )

ggsave('rank_urine_winners_in_stock_dot.pdf',
       path = '../plots/urine_tissue_manuscript')

### Calculate median or mean rank of all 40 top barcodes (10 per animal) in stock
urine_bc_levels %>%
    top_barcodes_by_max(10) %>%
    distinct(barcode) %>%
    left_join(cutoff_99pct_stock_barcodes %>%
              mutate(stock_rank = min_rank(desc(count)))) %>%
    pull(stock_rank) %>%
    median()

### Tables of median rank of urine winners in stock
urine_bc_levels %>%
    top_barcodes_by_max(10) %>%
    distinct(animal, barcode) %>%
    left_join(cutoff_99pct_stock_barcodes %>%
              mutate(stock_rank = min_rank(desc(count)))) %>%
    split(.$animal)

## -- Diversity Of Genomes In Urine ---------------------------------------

complexity_colors = c('Entropy' = 'cornflowerblue',
                      'Number To 50%' = 'navajowhite3',
                      'Number To 75%' = 'darkgoldenrod3',
                      'Diversity' = 'cornflowerblue')

urine_complexity_table = urine_bc_levels %>%
    group_by(animal, Days_pi) %>%
    arrange(desc(frac_weighted_count), .by_group=TRUE) %>%
    mutate(cum_frac = cumsum(frac_weighted_count)) %>%
    summarize(entropy = entropy(frac_weighted_count),
          num_to_50pct = sum(cum_frac < 0.5) + 1,
          num_to_75pct = sum(cum_frac < 0.75) + 1)

# Facet plot, but in log. In this case diversity is entropy
urine_complexity_table %>%
ggplot(aes(x=Days_pi)) +
         geom_line(aes(y = log(num_to_75pct), color = 'Number To 75%')) +
geom_line(aes(y = entropy, color = 'Entropy'), linewidth = 1) +
facet_wrap(vars(animal), ncol = 1, scales = 'free') +
         labs(title = bquote('Barcode Diversity In Urine'),
              subtitle = 'Log Scaled',
              x = 'Days Post-Injection') +
              ## y = 'Total Barcode Level') +
         scale_color_manual(name = 'Measure', values = complexity_colors[c('Entropy', 'Number To 75%')]) +
theme_minimal() +
theme(axis.text.y = element_blank(),
      axis.title.y = element_blank())

ggsave('diversity_entropy_num75_log.pdf',
               path='../plots/urine_tissue_manuscript')


## -- Urine Diversity Donuts ----------------------------------------------

plot_diversity_donut_w_colored_barcodes <- function(bc_table,
                                                    barcode_color_highlights,
                                                    faceting_variable) {

    winner_and_bg_colors = set_names(c(scales::hue_pal()(length(barcode_color_highlights)), 'gray75', 'gray60'),
                                     nm = c(barcode_color_highlights, 'background1', 'background2'))

    bc_table %>%
        group_by(sample) %>%
        arrange(desc(frac_weighted_count), .by_group=TRUE) %>%
        mutate(ymax = cumsum(frac_weighted_count)) %>%
        mutate(ymin = c(0, head(ymax, n=-1))) %>%
        ungroup() %>%
        select(animal, {{faceting_variable}}, barcode, frac_weighted_count, ymax, ymin) %>%
        mutate(bc_color = ifelse(row_number() %% 2, 'background1', 'background2')) %>%
        rowwise() %>%
        mutate(bc_color = ifelse(barcode %in% barcode_color_highlights,
                                barcode,
                                bc_color)) %>%
        ungroup() %>%
        arrange( {{faceting_variable}} ) %>%
        ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3.2, fill=bc_color)) +
        scale_fill_manual(values=winner_and_bg_colors, labels=names(winner_and_bg_colors)) +
        geom_rect(show.legend=FALSE) +
        coord_polar(theta='y') +
        xlim(c(2, 4)) +
        facet_wrap(vars( {{faceting_variable}} )) +
        ## labs(title = 'Working DL Diversity Donuts',
        ##      subtitle = 'Is this working?') +
        theme_classic() +
        theme(
            line = element_blank(),
            axis.line.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            ) +
        guides(fill = guide_legend('Barcode'))

    return(last_plot())
    }

winner_barcodes_by_animal = urine_bc_levels %>%
    top_barcodes_by_max(10) %>%
    group_by(animal) %>%
    distinct(barcode) %>%
    split(.$animal) %>%
    imap(~ pull(., barcode))

urine_bc_levels %>%
    split(.$animal) %>%
    imap(~ plot_diversity_donut_w_colored_barcodes(.x,
                                                   winner_barcodes_by_animal[[.y]],
                                                   Days_pi)) %>%
    imap(~ (.x + labs(title = 'Urine With Diversity Winners',
                      subtitle = paste('Mouse', .y)))) %>%
    iwalk(~ ggsave(paste0('urine_diversity_', .y, '_donut.pdf'),
                   plot = .x,
                   path = '../plots/urine_tissue_manuscript'))


## -- Stacked Area Charts For Winners--------------------------------------

plot_stacked_area_w_colored_barcodes <- function(bc_table, barcode_color_highlights) {

    winner_and_bg_colors = set_names(c(scales::hue_pal()(length(barcode_color_highlights)), 'gray75', 'gray60'),
                                     nm = c(barcode_color_highlights, 'background1', 'background2'))

    bc_table %>%
        ## filter(bc_level > 0) %>%
        group_by(animal, Days_pi) %>%
        ## arrange(desc(frac_weighted_count), .by_group=TRUE) %>%
        mutate(ymax = cumsum(frac_weighted_count)) %>%
        mutate(ymin = c(0, head(ymax, n=-1))) %>%
        ungroup() %>%
        select(animal, Days_pi, barcode, frac_weighted_count, ymax, ymin) %>%
        mutate(bc_color = ifelse(row_number() %% 2, 'background1', 'background2')) %>%
        rowwise() %>%
        mutate(bc_color = ifelse(barcode %in% barcode_color_highlights,
                                barcode,
                                bc_color)) %>%
        ungroup() %>%
        ## arrange(Days_pi) %>%
        ggplot(aes(x = Days_pi, ymax=ymax, ymin=ymin, fill=bc_color)) +
        scale_fill_manual(values=winner_and_bg_colors, labels=names(winner_and_bg_colors)) +
        geom_ribbon(show.legend=FALSE) +
        ## coord_polar(theta='y') +
        ## xlim(c(2, 4)) +
        ## labs(title = 'Working DL Diversity Donuts',
        ##      subtitle = 'Is this working?') +
        theme_classic() +
        theme(
            line = element_blank(),
            axis.line.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            ) +
        guides(fill = guide_legend('Barcode'))

    return(last_plot())
    }

winner_barcodes_by_animal = urine_bc_levels %>%
    top_barcodes_by_max(10) %>%
    group_by(animal) %>%
    distinct(barcode) %>%
    split(.$animal) %>%
    imap(~ pull(., barcode))

urine_bc_levels %>%
    filter(animal == 'FL') %>%
    plot_stacked_area_w_colored_barcodes(., winner_barcodes_by_animal)
    
ggsave('urine_winners_frac_FL_stacked_area.pdf',
       path = '../plots/urine_and_tissue')

urine_bc_levels %>%
    split(.$animal) %>%
    imap(~ plot_diversity_donut_w_colored_barcodes(.x, winner_barcodes_by_animal[[.y]])) %>%
    imap(~ (.x + labs(title = 'Urine With Diversity Winners',
                      subtitle = paste('Mouse', .y)))) %>%
    iwalk(~ ggsave(paste0('urine_diversity_', .y, '_donut.pdf'),
                   plot = .x,
                   path = '../plots/urine_and_tissue'))

## -- Relative Fold Increase For Shedding Events --------------------------

## If we're using median as baseline, excluding 0 from median calculation, make sure every barcode in every animal has at least one nonzero level. Zero medians would mess up ratio calculations
urine_bc_levels %>% group_by(animal, barcode) %>% summarize(max_bc_level = max(bc_level)) %>% filter(max_bc_level == 0)

median_bc_levels = urine_bc_levels %>%
    filter(bc_level > 0) %>%
    group_by(animal, barcode) %>%
    summarize(median_bc_level = median(bc_level))

max_ratio_vs_median_per_barcode = urine_bc_levels %>%
    left_join(median_bc_levels) %>%
    mutate(ratio_vs_median = bc_level / median_bc_level) %>%
    group_by(barcode) %>%
    slice_max(ratio_vs_median, n = 1) %>%
    ungroup() %>%
    mutate(log10_ratio_vs_median = log10(ratio_vs_median)) %>%
    select(-sample, -mean_ul_urine, -frac_weighted_count, -Days_pi) %>%
    arrange(desc(ratio_vs_median))

tibble(achievement = c('winner', 'slacker', 'loser'))

winner_slacker_loser_bc_levels = bind_rows(
    slice_barcodes_by_max(urine_bc_levels, 1, 10) %>%
    mutate(achievement = 'winner'),
    slice_barcodes_by_max(urine_bc_levels, 2001, 2010) %>%
    mutate(achievement = 'slacker'),
    slice_barcodes_by_max(urine_bc_levels, 4003, 4012) %>%
    mutate(achievement = 'loser')
)

winner_slacker_loser_median_ratios = winner_slacker_loser_bc_levels %>%
    left_join(median_bc_levels) %>%
    mutate(ratio_vs_median = bc_level / median_bc_level) %>%
    select(animal, Days_pi, barcode, bc_level, median_bc_level, ratio_vs_median, achievement)

## barcodes_grouped_by_achievement_fct = winner_slacker_loser_median_ratios %>%
##     filter(animal == 'FL') %>%
##     group_by(achievement) %>%
##     arrange(desc(ratio_vs_median)) %>%
##     distinct(barcode) %>%
##     pull(barcode) %>%
##     as_factor()

winner_slacker_loser_median_ratios %>%
    filter(animal == 'FL') %>%
    group_by(achievement) %>%
    arrange(desc(ratio_vs_median)) %>%
    ggplot(aes(x = Days_pi, y = fct_rev(as_factor(barcode)), height = log10(ratio_vs_median + 1), fill = achievement)) +
    geom_ridgeline(color = 'gray90',
                   scale = 0.5) +
    scale_fill_manual(name = 'Cohort', labels = c('Bottom 10', 'Middle 10', 'Top 10'),
                      values = c('goldenrod', 'darkolivegreen3', 'cornflowerblue'),
                      guide = guide_legend(reverse = TRUE)) +
    labs(title = '10-Fold Change vs Median Level Per Barcode',
         subtitle = 'Mouse FL',
         x = 'Days Post-Injection',
         y = 'Barcode') +
    theme_ridges()

ggsave('winner_slacker_loser_FL_ridge.pdf',
       path = '../plots/urine_and_tissue')

urine_bc_levels %>%
    top_barcodes_by_max(30) %>%
    split(.$animal) %>%
    imap(~ ggplot(.x, aes(x = Days_pi, y = barcode, height = bc_level)) +
             geom_ridgeline(fill = animal_colors[[.y]],
                            color = 'grey90',
                            scale = 0.0005) +
             scale_x_continuous(expand=c(0.01, 0)) +
             scale_y_discrete(expand=c(0.01, 0)) +
             labs(x = 'Days post-injection',
                  y = 'Barcode',
                  title = 'Top 30 Barcodes Ridge Plot',
                  subtitle = paste('Mouse', .y)) +
             theme_ridges()) %>%
    iwalk(~ ggsave(paste0('urine_top30_linear_', .y, '_ridge.pdf'),
                   plot = .x,
                   path = '../plots/urine_and_tissue',
                   height = 17))
             
ranked_urine_bc_levels = urine_bc_levels %>%
    group_by(animal, Days_pi) %>%
    arrange(desc(bc_level)) %>%
    mutate(rank = 1:n())

winner_slacker_loser_ranks = winner_slacker_loser_bc_levels %>%
    left_join(ranked_urine_bc_levels)

winner_slacker_loser_ranks %>%
    filter(animal == 'FL') %>%
    group_by(achievement) %>%
    ## arrange(rank) %>%
    ggplot(aes(x = Days_pi, y = fct_rev(as_factor(barcode)), height = 4013 - rank, fill = achievement)) +
    geom_ridgeline(color = 'gray90',
                   scale = 0.0002) +
    scale_fill_manual(name = 'Cohort', labels = c('Bottom 10', 'Middle 10', 'Top 10'),
                      values = c('goldenrod', 'darkolivegreen3', 'cornflowerblue'),
                      guide = guide_legend(reverse = TRUE)) +
    labs(title = 'Rank Of Barcodes Per Sample',
         subtitle = 'Mouse FL',
         x = 'Days Post-Injection',
         y = 'Barcode') +
    theme_ridges()

ggsave('ranks_winner_slacker_loser_FL_ridge.pdf',
       path = '../plots/urine_and_tissue')


## -- Stock vs Urine Composition UpSet Plots-------------------------------

urine_bc_levels %>%
    filter(frac_weighted_count >= 0.0001) %>%
    group_by(animal, Days_pi) %>%
    summarize(n = n()) %>%
    ggplot(aes(x = Days_pi, y = n)) +
    geom_line(aes(color = as_factor(animal)), show.legend = FALSE) +
    facet_wrap(vars(animal), ncol = 1) +
    labs(title = 'Number Of Barcodes With More Than 0.01% Expression',
         x = 'Days Post-Injection',
         y = 'Number of Barcodes') +
    scale_color_manual(values = animal_colors) +
    theme_minimal()
    
ggsave('num_bcs_gt_0.01pct_line.pdf',
       path = '../plots/urine_and_tissue')

## stock_urine_lte20_gt20 = list('Stock' = (cutoff_99pct_stock_barcodes %>% pull(barcode)),
##                            'Urine Up To Day 20' = (urine_bc_levels %>% filter(Days_pi <= 20) %>% pull(barcode)),
##                            'Urine After Day 20' = (urine_bc_levels %>% filter(Days_pi > 20) %>% pull(barcode)))

stock_urine_lte20_gt20 = list('Stock' = cutoff_99pct_stock_barcodes,
                           'Urine Up To Day 20' = filter(urine_bc_levels, Days_pi <= 20 & bc_level > 0),
                           'Urine After Day 20' = filter(urine_bc_levels, Days_pi > 20 & bc_level > 0)) %>%
    map(~ pull(., barcode)) %>%
    map(~ unique(.))


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
    mutate(library = if_else(library == 'Ligated', 'Ligated virus genomes', paste(library, 'library'))) %>%
    mutate(intersection = 1) %>%
    pivot_wider(id_cols = c(library, barcode, viral_count),
                names_from = library,
                values_from = intersection) %>%
    replace_na(list(`Plasmid library`= 0, `Ligated virus genomes` = 0, `Viral library` = 0)) %>%
    upset(c('Viral library', 'Ligated virus genomes', 'Plasmid library'),
          name = 'Intersection Of Clustered Barcodes',
          set_sizes = FALSE,
          sort_sets = FALSE,
          sort_intersections = 'ascending',
          sort_intersections_by = 'degree',
          intersections = list(
              c('Plasmid library', 'Viral library'),
              c('Ligated virus genomes', 'Viral library'),
              c('Plasmid library', 'Ligated virus genomes', 'Viral library'),
              'Viral library'
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


## -- Decline Of Distinct Barcodes Over Time ------------------------------

urine_bc_levels %>%
    group_by(animal, barcode) %>%
    filter(bc_level > 0) %>%
    summarize(last_day = max(Days_pi)) %>%
    count(last_day) %>% filter(animal == 'MR')


# --- Tissue Analysis -----------------------------------------------------

tissue_pcr = read_tsv('organs_qPCR_table.tsv', na = c('N/A')) %>%
    rename(sample=Sample_Name, animal=Animal, organ=Organ,
           quantity_mean_ug_dna=`Quantity_Mean/ug_DNA`,
           quantity_mean_ul_blood=`Quantity_Mean/ul_blood_or_virus_stock`) %>%
    select(sample, animal, organ, quantity_mean_ug_dna, quantity_mean_ul_blood) %>%
    mutate(blood = str_detect(organ, '.*blood.*')) %>%
    mutate(organ = recode(organ,
                          Gonades = 'Gonads',
                          `Plasma*` = 'Plasma',
                          Testicules = 'Testicles',
                          `Whole blood*` = 'Whole blood'))

tissue_bc_levels = tissue_to_stock_barcodes %>%
    group_by(sample) %>%
    mutate(frac_weighted_count = weighted_count / sum(weighted_count)) %>%
    ungroup() %>%
    left_join(tissue_pcr) %>%
    mutate(bc_level = if_else(blood,
                              frac_weighted_count * quantity_mean_ul_blood * 1000,
                              frac_weighted_count * quantity_mean_ug_dna)) %>%
    select(sample, animal, organ, stock_bc, bc_level, frac_weighted_count, quantity_mean_ug_dna, quantity_mean_ul_blood, blood) %>%
    rename(barcode = stock_bc) %>%
    filter(sample != 'PLMR')        # bad sample

## -- Rank Tissue Winners In Stock ----------------------------------------

tissue_bc_levels %>%
    top_barcodes_by_max(10) %>%
    filter(bc_level == max_level) %>%
    split(.$animal) %>%
    imap(~ mutate(.x, tissue_rank = desc(max_level) %>% min_rank() %>% as_factor())) %>%
    bind_rows() %>% 
    left_join(cutoff_99pct_stock_barcodes %>%
              mutate(stock_rank = desc(count) %>% min_rank())) %>%
    select(animal, max_level, barcode, tissue_rank, count, stock_rank) %>%
    ggplot(aes(x = fct_relevel(tissue_rank, rev), y = stock_rank)) + 
    geom_point(aes(color = as_factor(animal)), size=3, show.legend=FALSE) +
    facet_wrap(vars(animal)) +
    labs(title = 'Rank Of Tissue Winners In Stock',
         x = "Rank Of Winner In Tissue",
         y = "Rank In Stock") +
    scale_y_continuous(limits = c(1, dim(cutoff_99pct_stock_barcodes)[[1]])) +
    scale_x_discrete(breaks = seq(1, 10)) +
    scale_color_manual(values=animal_colors) +
    coord_flip() +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(linetype = 'dotted',
                                            color = 'gray70'),
          panel.grid.minor.y = element_blank(),
          panel.spacing = unit(2, 'lines'),
          strip.background = element_rect(fill = 'gray90',
                                          color = 'gray90')
          )

ggsave('rank_tissue_winners_in_stock_dot.pdf',
       path = '../plots/urine_tissue_manuscript')


## MWU to test if tissue winners are overrepresented in stock

bladder_kidney_winners = tissue_bc_levels %>%
    filter(organ %in% c('Bladder', 'Kidney')) %>%
    top_barcodes_by_max(10) %>%
    pull(barcode) %>%
    unique()

bladder_kidney_winners_within_stock = cutoff_99pct_stock_barcodes %>%
    mutate(is_winner = if_else(barcode %in% bladder_kidney_winners, 'winner', 'not_winner'))

wilcox.test(ranking ~ is_winner, data = bladder_winners_within_stock,


## -- Correlation Of Total Shed For Tissue vs Urine -----------------------

cor(urine_bc_levels %>%
    filter(animal == 'FL') %>%
    top_barcodes_by_max(100) %>%
    group_by(barcode) %>%
    summarize(total_shed = sum(bc_level)) %>%
    slice_max(total_shed, n = 400) %>%
    arrange(barcode) %>%
    pull(total_shed),
    tissue_bc_levels %>%
    filter(animal == 'FL') %>%
    top_barcodes_by_max(100) %>%
    group_by(barcode) %>%
    summarize(total_shed = sum(bc_level)) %>%
    slice_max(total_shed, n = 400) %>%
    arrange(barcode) %>% # Want to compare the same barcode in both sets
    pull(total_shed), 
    method = 'spearman')  

cor(urine_bc_levels %>%
    filter(animal == 'FL') %>%
    top_barcodes_by_max(100) %>%
    group_by(barcode) %>%
    top_n(1, wt = bc_level) %>%
    

tissue_FL_10_winners = tissue_bc_levels %>%
    filter(animal == 'FL') %>%
    top_barcodes_by_max(10) %>%
    distinct(barcode) %>%
    pull(barcode)

cor(urine_bc_levels %>%
    filter(animal == 'FL' & barcode %in% tissue_FL_10_winners) %>%
    group_by(barcode) %>%
    summarize(total_shed = sum(bc_level)) %>%
    pull(total_shed),
    tissue_bc_levels %>%
    filter(animal == 'FL' & barcode %in% tissue_FL_10_winners) %>%
    group_by(barcode) %>%
    summarize(total_shed = sum(bc_level)) %>%
    pull(total_shed),
    method = 'spearman')

# Try taking 2-5 from each organ, rather than just top 10 overall?

tissue_FL_5_from_each_organ = tissue_bc_levels %>%
    filter(animal == 'FL') %>%
    group_by(organ) %>%
    slice_max(bc_level, n = 5) 


tissue_FL_10_from_each_organ = tissue_bc_levels %>%
    filter(animal == 'FL') %>%
    group_by(organ) %>%
    slice_max(bc_level, n = 10) %>%
    mutate(organ_rank = desc(bc_level) %>% min_rank() %>% as_factor()) %>%
    select(-quantity_mean_ul_blood, -frac_weighted_count)

tissue_FL_10_from_each_organ %>%
    left_join(urine_bc_levels %>%
            filter(animal == 'FL') %>%
            group_by(barcode) %>%
            summarize(total_shed = sum(bc_level)) %>%
            mutate(urine_rank = desc(total_shed) %>% min_rank())) %>%
    select(organ, barcode, urine_rank, organ_rank) %>%
    ggplot(aes(x = urine_rank, y = fct_relevel(organ_rank, rev))) +
    geom_point(color = animal_colors[['FL']], size = 3, show.legend = FALSE) +
    facet_wrap(vars(organ)) +
    scale_x_continuous(limits = c(1, nrow(distinct(urine_bc_levels, barcode)))) +
    labs(title = 'Rank Of Organ Winners In Urine',
         x = "Rank Of Winner In Urine",
         y = "Rank Of Winner In Organ") +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(linetype = 'dotted',
                                            color = 'gray70'),
          panel.grid.minor.y = element_blank(),
          panel.spacing = unit(2, 'lines'),
          strip.background = element_rect(fill = 'gray90',
                                          color = 'gray90')
          )

ggsave('rank_FL_organ_winners_in_urine_by_organ_dot.pdf',
       path = '../plots/urine_tissue_manuscript',
       height = 11,
       width = 10)

## Now for all animals

tissue_10_from_each_organ_by_animal = tissue_bc_levels %>%
    split(.$animal) %>%
    imap(~ group_by(.x, organ) %>%
         slice_max(bc_level, n = 10) %>%
         mutate(organ_rank = desc(bc_level) %>% min_rank() %>% as_factor()) %>%
         select(-quantity_mean_ul_blood, -frac_weighted_count)) 

urine_total_shed_by_animal = urine_bc_levels %>%
    split(.$animal) %>%
    imap(~ group_by(.x, barcode) %>%
             summarize(total_shed = sum(bc_level)) %>%
             mutate(urine_rank = desc(total_shed) %>% min_rank()))

map2_dfr(tissue_10_from_each_organ_by_animal, urine_total_shed_by_animal, left_join) %>%
split(.$animal) %>%
    imap(~ ggplot(.x, aes(x = urine_rank, y = fct_relevel(organ_rank, rev))) +
    geom_point(color = animal_colors[[.y]], size = 3, show.legend = FALSE) +
    facet_wrap(vars(organ)) +
    scale_x_continuous(limits = c(1, nrow(distinct(urine_bc_levels, barcode)))) +
    labs(title = paste0('Rank Of Organ Winners In Urine For ', .y),
         x = "Rank Of Winner In Urine",
         y = "Rank Of Winner In Organ") +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(linetype = 'dotted',
                                            color = 'gray70'),
          panel.grid.minor.y = element_blank(),
          panel.spacing = unit(2, 'lines'),
          strip.background = element_rect(fill = 'gray90',
                                          color = 'gray90')
          )
    ) %>%
iwalk(~ ggsave(paste('rank', .y, 'organ_winners_in_urine_by_organ_dot.pdf', sep = '_'),
               plot = .x,
               path = '../plots/urine_tissue_manuscript',
               height = 11,
               width = 10))

## Rank in last day urine, rather than total shed

tissue_10_from_each_organ_by_animal = tissue_bc_levels %>%
    split(.$animal) %>%
    imap(~ group_by(.x, organ) %>%
         slice_max(bc_level, n = 10) %>%
         mutate(organ_rank = desc(bc_level) %>% min_rank() %>% as_factor()) %>%
         select(-quantity_mean_ul_blood, -frac_weighted_count)) 

urine_last_day_rank_by_animal = urine_bc_levels %>%
    group_by(animal) %>%
    mutate(last_day = max(Days_pi)) %>%
    filter(Days_pi == last_day) %>% 
    ungroup() %>%
    split(.$animal) %>%
    imap(~ mutate(.x, last_day_rank = desc(bc_level) %>% min_rank())) %>%
    imap(~ select(.x, barcode, last_day_rank))

map2_dfr(tissue_10_from_each_organ_by_animal, urine_last_day_rank_by_animal, left_join) %>%
    split(.$animal) %>%
    imap(~ ggplot(.x, aes(x = last_day_rank, y = fct_relevel(organ_rank, rev))) +
    geom_point(color = animal_colors[[.y]], size = 3, show.legend = FALSE) +
    facet_wrap(vars(organ)) +
    scale_x_continuous(limits = c(1, nrow(distinct(urine_bc_levels, barcode)))) +
    labs(title = paste0('Rank Of Organ Winners In Last Day Of Urine For ', .y),
         x = "Rank Of Winner In Urine",
         y = "Rank Of Winner In Organ") +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(linetype = 'dotted',
                                            color = 'gray70'),
          panel.grid.minor.y = element_blank(),
          panel.spacing = unit(2, 'lines'),
          strip.background = element_rect(fill = 'gray90',
                                          color = 'gray90')
          )
    ) %>%
iwalk(~ ggsave(paste('rank', .y, 'organ_winners_in_last_day_urine_by_organ_dot.pdf', sep = '_'),
               plot = .x,
               path = '../plots/urine_tissue_manuscript',
               height = 11,
               width = 10))

## -- Are Urine Winners In Bottom Quartile Of Stock Overrepresented In Tissue ?

cutoff_99pct_stock_barcodes %>%
    pull(count) %>%
    quantile(., prob = c(0.25, 0.5, 0.75))

#     25%     50%     75%
# 3184.75 4944.00 8037.75

bottom_quartile_stock_barcodes = cutoff_99pct_stock_barcodes %>%
    filter(count < 3184.75)

bottom_quartile_stock_cutoff = cutoff_99pct_stock_barcodes %>%
    pull(count) %>%
    quantile(., prob = c(0.25, 0.5, 0.75)) %>%
    pluck('25%')

bottom_quartile_stock_barcodes = cutoff_99pct_stock_barcodes %>%
    filter(count < bottom_quartile_stock_cutoff)

# Let's look at top 5% of tissue barcodes (this is across all mice)

tissue_bc_levels %>%
    pull(bc_level) %>%
    quantile(., prob = c(0.95, 1))

#          95%         100%
# 1.020217e+02 2.352305e+08

top_5pct_tissue_bcs_all_mice = tissue_bc_levels %>%
    filter(bc_level > 1.020217e+02) %>%
    distinct(barcode)

# Let's look at top 5% of urine barcodes (this is across all mice)

urine_bc_levels %>%
    pull(bc_level) %>%
    quantile(., prob = c(0.95, 1))

#          95%         100%
# 3.169872e-02 1.495440e+04

top_5pct_urine_bcs_all_mice = urine_bc_levels %>%
    filter(bc_level > 3.169872e-02) %>%
    distinct(barcode)

top_5pct_urine_cutoffs = urine_bc_levels %>%
    split(.$animal) %>%
    map(~ group_by(., barcode)) %>%  # bc_level has a value for every time point, choose just the one with highest value
    map(~ slice_max(., order_by = bc_level)) %>%
    map(~ ungroup(.)) %>%
    map(~ pull(., bc_level)) %>%
    map(~ quantile(., prob = c(0.95, 1))) %>%
    map(~ pluck(., '95%'))

top_5pct_urine_all_animals = urine_bc_levels %>%
    split(.$animal) %>%
    map(~ group_by(., barcode)) %>%
    map(~ slice_max(., order_by = bc_level)) %>%
    map(~ ungroup(.)) %>%
    map2(., top_5pct_urine_cutoffs, ~ filter(.x, bc_level > .y)) 

top_5pct_urine_in_bottom_stock = top_5pct_urine_all_animals %>%
    map(~ filter(., barcode %in% pull(bottom_quartile_stock_barcodes, barcode)))

top_5pct_tissue_cutoffs = tissue_bc_levels %>%
    split(.$animal) %>%
    map(~ group_by(., barcode)) %>%  # bc_level has a value for every time point, choose just the one with highest value
    map(~ slice_max(., order_by = bc_level)) %>%
    map(~ ungroup(.)) %>%
    map(~ pull(., bc_level)) %>%
    map(~ quantile(., prob = c(0.95, 1))) %>%
    map(~ pluck(., '95%'))

top_5pct_tissue_all_animals = tissue_bc_levels %>%
    split(.$animal) %>%
    map(~ group_by(., barcode)) %>%
    map(~ slice_max(., order_by = bc_level)) %>%
    map(~ ungroup(.)) %>%
    map2(., top_5pct_tissue_cutoffs, ~ filter(.x, bc_level > .y)) 

top_5pct_tissue_in_bottom_stock = top_5pct_tissue_all_animals %>%
    map(~ filter(., barcode %in% pull(bottom_quartile_stock_barcodes, barcode)))

map2(top_5pct_urine_all_animals, top_5pct_tissue_all_animals, ~ intersect(pull(.x, barcode), pull(.y, barcode))) %>%
    map(~ length(.))

# 2nd attempt

bottom_quartile_stock_barcodes = cutoff_99pct_stock_barcodes %>%
    top_frac(-0.25, wt = count)

max_urine_bc_levels = urine_bc_levels %>%
    group_by(animal, barcode) %>%
    summarize(max_level = max(bc_level)) %>%
    mutate(total_level = sum(max_level), frac = max_level / total_level) %>%
    top_frac(0.05, wt = max_level)

top_5pct_urine_all_animals = urine_bc_levels %>%
    group_by(animal, barcode) %>%
    slice_max(order_by = bc_level) %>%
    ungroup(barcode) %>%
    top_frac(0.05, wt = bc_level) %>%
    ungroup()

top_5pct_tissue = tissue_bc_levels %>%
    group_by(animal, organ) %>%
    top_frac(0.05, wt = bc_level) %>%
    ungroup(animal)

inner_join(top_5pct_tissue, top_5pct_urine_all_animals, by = c('animal', 'barcode'))

inner_join(top_5pct_tissue, top_5pct_urine_all_animals, by = c('animal', 'barcode')) %>% count(animal, organ) %>% slice_max(order_by = n)

max_tissue_frac_weighted_count = tissue_bc_levels %>%
    group_by(animal, barcode) %>%
    
top_5pct_tissue = tissue_bc_levels %>%
    group_by(animal, barcode) %>%
    slice_max(order_by = frac_weighted_count) %>%
    ungroup(barcode) %>%
    top_frac(0.05, wt = frac_weighted_count) %>%
    ungroup()

fisher_top5pct_pvalues = c()
for (animal_name in names(animal_colors)) {
    top_5pct_urine_top_5pct_tissue = intersect(
        filter(top_5pct_urine_all_animals, animal == animal_name) %>% pull(barcode),
        filter(top_5pct_tissue, animal == animal_name) %>% pull(barcode)) %>%
        length()
    top_5pct_urine_bot_95pct_tissue = 200 - top_5pct_urine_top_5pct_tissue
    barcode_top_5pct_quintile_bins = rbind(c(top_5pct_urine_top_5pct_tissue, top_5pct_urine_bot_95pct_tissue),
                                           c(top_5pct_urine_bot_95pct_tissue, 4012 - 2 * top_5pct_urine_bot_95pct_tissue - top_5pct_urine_top_5pct_tissue))
    fisher_top5pct_pvalues = c(fisher_top5pct_pvalues, fisher.test(barcode_top_5pct_quintile_bins)$p.value)
    }
names(fisher_top5pct_pvalues) = names(animal_colors)


## -- Urine Winners In Stock And Tissue -----------------------------------

## Get median rank of urine in stock
urine_rank_in_stock_table = urine_bc_levels %>%
    split(.$animal) %>%
    imap(~ top_barcodes_by_max(.x, 10)) %>%
    imap(~ arrange(.x, bc_level)) %>%
    imap(~ left_join(.x, mutate(cutoff_99pct_stock_barcodes, stock_ranking = desc(count) %>% min_rank()))) %>%
    bind_rows() %>%
    group_by(animal) %>%
    summarize(median_stock_rank = median(stock_ranking)) %>%
    rename(Animal = animal, Urine = median_stock_rank)

## %>%
##     rename(Animal = animal, `Median Rank` = median_stock_rank) %>%
##     kbl(caption = 'Median Rank Of Urine Winners In Stock') %>%
##     kable_styling(bootstrap_options='striped', full_width=FALSE) %>%
##     save_kable('../plots/urine_tissue_manuscript/urine_winner_median_stock_ranks_table.html')

## Get median rank of tissue in stock
tissue_rank_in_stock_table = tissue_bc_levels %>%
    split(.$animal) %>%
    imap(~ top_barcodes_by_max(.x, 10)) %>%
    imap(~ arrange(.x, bc_level)) %>%
    imap(~ left_join(.x, mutate(cutoff_99pct_stock_barcodes, stock_ranking = desc(count) %>% min_rank()))) %>%
    bind_rows() %>%
    group_by(animal) %>%
    summarize(median_stock_rank = median(stock_ranking)) %>%
    rename(Animal = animal, Tissue = median_stock_rank)

## %>%
##     rename(Animal = animal, `Median Rank` = median_stock_rank) %>%
##     kbl(caption = 'Median Rank Of Tissue Winners In Stock') %>%
##     kable_styling(bootstrap_options='striped', full_width=FALSE) %>%
##     save_kable('../plots/urine_tissue_manuscript/tissue_winner_median_stock_ranks_table.html')

left_join(urine_rank_in_stock_table, tissue_rank_in_stock_table) %>%
    kbl(caption = 'Median Rank Of Urine And Tissue Winners In Stock') %>%
    kable_styling(bootstrap_options = 'striped', full_width = FALSE) %>%
    column_spec(2, width = '8em') %>%
    column_spec(3, width = '8em') %>%
    save_kable('../plots/urine_tissue_manuscript/urine_and_tissue_winner_median_stock_ranks_table.html')

## Overlap of tissue and urine winners

each_animal_urine_winners = urine_bc_levels %>% split(.$animal) %>% imap(~ top_barcodes_by_max(.x, 10)) %>% imap(~ distinct(.x, barcode)) %>% imap(~ pull(.x, barcode))

each_animal_tissue_winners = tissue_bc_levels %>% split(.$animal) %>% imap(~ top_barcodes_by_max(.x, 10)) %>% imap(~ distinct(.x, barcode)) %>% imap(~ pull(.x, barcode))

map2(each_animal_urine_winners, each_animal_tissue_winners, \(x, y) intersect(x, y))

## -- Tissue Winners In Urine ---------------------------------------------

tissue_FL_10_winners = tissue_bc_levels %>%
    filter(animal == 'FL') %>%
    top_barcodes_by_max(10) %>%
    distinct(barcode) %>%
    pull(barcode)

tissue_FL_50_winners = tissue_bc_levels %>%
    filter(animal == 'FL') %>%
    top_barcodes_by_max(50) %>%
    distinct(barcode) %>%
    pull(barcode)

urine_bc_levels %>%
    filter(animal == 'FL' & barcode %in% tissue_FL_50_winners) %>%
    arrange(desc(bc_level))

urine_bc_levels %>%
    filter(animal == 'FL' & barcode %in% tissue_FL_10_winners) %>%
    ggplot(aes(x = Days_pi, fill = barcode)) +
    geom_area(aes(y = mean_ul_urine), fill = 'grey85') +
    geom_area(aes(y = bc_level)) +
                facet_wrap(vars(barcode), ncol = 1) +
    theme_minimal() +
                theme(legend.position = 'none')

ggsave('tissue_winner_urine_plots_FL_area.pdf',
       path = '../plots/urine_tissue_manuscript',
       height = 10)

## -- GC Content Top 10 Tissues Per Mouse ---------------------------------

gc_pct_all_tissue_barcodes = tissue_bc_levels %>%
    distinct(barcode) %>%
    mutate(GC_pct = str_count(barcode, '[GC]') / str_length(barcode)) %>%
    add_column(animal = 'All Barcodes', .before = 1)

tissue_bc_levels %>%
    top_barcodes_by_max(10) %>%
    distinct(animal, barcode) %>%
    mutate(GC_pct = str_count(barcode, '[GC]') / str_length(barcode)) %>%
    bind_rows(., gc_pct_all_barcodes) %>%
    mutate(animal = fct_relevel(animal, c('FL', 'FR', 'ML', 'MR', 'All Barcodes'))) %>%
    ggplot(aes(animal, GC_pct)) +
    geom_boxplot(aes(color = animal), linewidth = 0.8, show.legend = FALSE) +
    labs(title = 'GC Content For Top 10 Barcodes Per Animal In Urine',
         x = 'Mouse',
         y = 'GC Content') +
    scale_color_manual(values = c(animal_colors, `All Barcodes` = 'gray70')) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format()) +
    theme_minimal()

ggsave('GC_urine_winners_boxplot.pdf',
        path = '../plots/urine_tissue_manuscript')

tissue_bc_levels %>%
    slice_max(n = 10, order_by = bc_level, by = c(animal, organ)) %>%
    group_by(animal) %>%
    distinct(barcode) %>%
    mutate(GC_pct = str_count(barcode, '[GC]') / str_length(barcode)) %>%
    ungroup() %>%
    bind_rows(., gc_pct_all_barcodes) %>%
    mutate(animal = fct_relevel(animal, c('FL', 'FR', 'ML', 'MR', 'All Barcodes'))) %>%
    ggplot(aes(animal, GC_pct)) +
    geom_boxplot(aes(color = animal), linewidth = 0.8, show.legend = FALSE) +
    labs(title = 'GC Content For Top 10 Barcodes Per Animal In Urine',
         x = 'Mouse',
         y = 'GC Content') +
    scale_color_manual(values = c(animal_colors, `All Barcodes` = 'gray70')) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format()) +
    theme_minimal()

ggsave('GC_tissue_winners_boxplot.pdf',
        path = '../plots/urine_tissue_manuscript')

## -- Tissue Diversity Donuts ---------------------------------------------

winner_barcodes_by_animal = urine_bc_levels %>%
    top_barcodes_by_max(10) %>%
    group_by(animal) %>%
    distinct(barcode) %>%
    split(.$animal) %>%
    imap(~ pull(., barcode))

tissue_bc_levels %>%
    split(.$animal) %>%
    imap(~ plot_diversity_donut_w_colored_barcodes(.x,
                                                   winner_barcodes_by_animal[[.y]],
                                                   organ)) %>%
    imap(~ (.x + labs(title = 'Tissue Diversity With Urine Winners',
                      subtitle = paste('Mouse', .y)))) %>%
    iwalk(~ ggsave(paste0('tissue_diversity_w_urine_winners_', .y, '_donut.pdf'),
                   plot = .x,
                   path = '../plots/urine_and_tissue'))


## -- Tissue Heatmaps -----------------------------------------------------

save_pheatmap_pdf <- function(x, path, filename, width, height) {
    pdf(paste(path, filename, sep='/'), w=width, h=height, onefile=FALSE)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
    }

scaled_matrix_tissue_bc_levels_list = tissue_bc_levels %>%
    select(animal, organ, barcode, bc_level) %>%
    pivot_wider(names_from = organ, values_from = bc_level, values_fill = 0) %>%
    mutate_if(is.numeric, replace_na, 0) %>%
    split(.$animal) %>%
    imap(~ select(., -animal)) %>%
    imap(~ column_to_rownames(., var = 'barcode')) %>%
    imap(~ as.matrix(.x))

%>%
    imap(~ scale(.x))  # to z-score

selected_bcs_for_heatmaps = tissue_bc_levels %>%
    top_barcodes_by_max(30) %>%
    select(animal, barcode) %>%
    distinct() %>%
    group_by(animal) %>%
    nest() %>%
    mutate(data = map(data, ~ pull(.x, barcode))) %>%
    deframe()

scaled_matrix_tissue_bc_levels_list = tissue_bc_levels %>%
    split(.$animal) %>%
    imap(~ select(.x, barcode, organ, bc_level)) %>%
    imap(~ filter(.x, barcode %in% selected_bcs_for_heatmaps[[.y]])) %>%
    imap(~ pivot_wider(.x, names_from = organ, values_from = bc_level)) %>%
    imap(~ column_to_rownames(.x, var = 'barcode')) %>%
    imap(~ log10(.x + 0.001)) %>%
    imap(~ scale(.x)) 

## FL_selected_barcodes = tissue_bc_levels %>%
##     filter(animal == 'FL') %>%
##     group_by(barcode) %>%
##     summarize(variance = var(bc_level)) %>%
##     slice_max(order_by = variance, n = 30)

## FL_scaled_matrix = tissue_bc_levels %>%
##     filter(animal == 'FL' & barcode %in% pull(FL_selected_barcodes, barcode)) %>%
##     select(organ, barcode, bc_level) %>%
##     pivot_wider(names_from = organ, values_from = bc_level, values_fill = 0) %>%
##     column_to_rownames('barcode') %>%
##     as.matrix()

for (animal in c('FL', 'FR', 'ML', 'MR')) {
    FL_scaled_heatmap = pheatmap(scaled_matrix_tissue_bc_levels_list[[animal]],
                                 cluster_rows = FALSE,
                                 main=paste0('Tissue: ', animal, ' Top 30 Urine Winners z-score Of log10'))
    save_pheatmap_pdf(FL_scaled_heatmap, '../plots/urine_and_tissue', paste('bc_level_30', animal, 'scaled_heatmap.pdf', sep = '_'), height = 12, width = 12)
    }
