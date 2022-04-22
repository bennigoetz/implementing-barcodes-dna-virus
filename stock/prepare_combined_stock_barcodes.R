library(tidyverse)

clustered_combined_stock_barcodes = read_tsv('BC7_combined_stock_mp_L3_clusters.tsv',
                                             col_names = c('barcode', 'count', 'elements'),
                                             progress = FALSE) %>%
    select(-elements)

cutoff_99pct_stock_barcodes = clustered_combined_stock_barcodes %>%
    mutate(total = sum(count),
           pct = count / total,
           cum_pct = cumsum(pct)) %>%
    filter(cum_pct <= 0.99)
