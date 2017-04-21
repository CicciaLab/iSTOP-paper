library(tidyverse)
library(iSTOP)

# Change the file paths as necessary, the files are on our shared Google Drive
COSMIC <- read_csv('~/Google Drive/Rothstein-Lab/2017-iSTOP/Eric/data/COSMIC-nonsense.csv') %>% rename(genome_coord = coord)
CDS    <- read_csv('~/Google Drive/Rothstein-Lab/2017-iSTOP/Eric/data/CDS-Human.csv')

my_gene = 'BRCA1'

iSTOP <-
  filter(CDS, gene == my_gene) %>%
  locate_codons(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  locate_iSTOP(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  # Within each gene filter for longest isoforms
  group_by(gene) %>%
#  filter(cds_length == max(cds_length)) %>%
#  add_RFLP(width = 150) %>% # you can change the width as you please (150 is the default and the max allowed)
#  add_RFLP(width = 100) %>% # you can add as many RFLP widths as you like
  add_RFLP(width = 50)

isoforms <- plot_spliced_isoforms(
  gene   = my_gene,
  coords = filter(CDS, tx %in% iSTOP$tx),
  colors = c('red', 'black', 'blue', 'darkgreen'),
  `Nonsense in cancer`   = COSMIC,
  `CAA, CAG, CGA, TGG`   = iSTOP,
  `iSTOP targetable`     = filter(iSTOP, match_any),
  `Verifiable with RFLP` = filter(iSTOP, match_any & has(RFLP_50)) # The RFLP number will change depending on the width
)
isoforms
ggsave('~/Desktop/SMARCAL1-2865.pdf', isoforms, width = 8, height = 3)
