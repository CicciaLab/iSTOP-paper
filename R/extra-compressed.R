library(tidyverse)
library(stringr)

fix_NMD <- function(file) {
  read_csv(file, col_types = cols(percent_tx = col_double())) %>%
    mutate(
      NMD_pred_tmp = str_replace_all(NMD_pred, 'TRUE', 'false'),
      NMD_pred_tmp = str_replace_all(NMD_pred_tmp, 'FALSE', 'TRUE'),
      NMD_pred     = str_replace_all(NMD_pred_tmp, 'false', 'FALSE')
    ) %>%
    select(-NMD_pred_tmp) %>%
    write_csv(file.path('data/iSTOP-compact-new', basename(file)))
}


compress_further <- function(file) {
  read_csv(file, col_types = cols(percent_tx = col_double())) %>%
    select(-searched) %>%
    write_csv(file.path('data/iSTOP-compact-new', basename(file)))
}

list.files('data/iSTOP-compact', '[.]csv$', full.names = T) %>%
  walk(~compress_further(.))
