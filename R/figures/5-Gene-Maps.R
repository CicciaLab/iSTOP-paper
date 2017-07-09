library(tidyverse)
library(iSTOP)

rm(list = ls())

# Change the file paths as necessary, the files are on our shared Google Drive
COSMIC <- read_csv('data/COSMIC/COSMIC-nonsense.csv')
CDS    <- read_csv('data/CDS/Hsapiens-hg38.csv')

# ---- FANCM ----
iSTOP_FANCM <-
  filter(CDS, gene == 'FANCM') %>%
  locate_codons(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  locate_PAM(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  group_by(gene) %>%
  filter(cds_length == max(cds_length)) %>%
  add_RFLP(width = 50)

FANCM <- plot_spliced_isoforms(
  gene   = 'FANCM',
  coords = filter(CDS, tx %in% iSTOP_FANCM$tx),
  colors = c('black', 'blue', 'darkgreen'),
  `CAA, CAG, CGA, TGG`   = iSTOP_FANCM,
  `iSTOP targetable`     = filter(iSTOP_FANCM, match_any),
  `Verifiable with RFLP` = filter(iSTOP_FANCM, match_any & has(RFLP_50))
)
FANCM
ggsave('figures/FANCM.pdf', FANCM, width = 8, height = 3)


# ---- TIMELESS ----
iSTOP_TIMELESS <-
  filter(CDS, gene == 'TIMELESS') %>%
  locate_codons(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  locate_PAM(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  group_by(gene) %>%
  filter(cds_length == max(cds_length)) %>%
  add_RFLP(width = 50)

TIMELESS <- plot_spliced_isoforms(
  gene   = 'TIMELESS',
  coords = filter(CDS, tx %in% iSTOP_TIMELESS$tx),
  colors = c('black', 'blue', 'darkgreen'),
  #`Nonsense in cancer`   = COSMIC,
  `CAA, CAG, CGA, TGG`   = iSTOP_TIMELESS,
  `iSTOP targetable`     = filter(iSTOP_TIMELESS, match_any),
  `Verifiable with RFLP` = filter(iSTOP_TIMELESS, match_any & has(RFLP_50)) # The RFLP number will change depending on the width
)
TIMELESS
ggsave('figures/TIMELESS.pdf', TIMELESS, width = 8, height = 3)

# ---- SPRTN ----
iSTOP_SPRTN <-
  filter(CDS, gene == 'SPRTN') %>%
  locate_codons(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  locate_PAM(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  group_by(gene) %>%
  filter(cds_length == max(cds_length)) %>%
  add_RFLP(width = 50)

SPRTN <- plot_spliced_isoforms(
  gene   = 'SPRTN',
  coords = filter(CDS, tx %in% iSTOP_SPRTN$tx),
  colors = c('black', 'blue', 'darkgreen'),
  #`Nonsense in cancer`   = COSMIC,
  `CAA, CAG, CGA, TGG`   = iSTOP_SPRTN,
  `iSTOP targetable`     = filter(iSTOP_SPRTN, match_any),
  `Verifiable with RFLP` = filter(iSTOP_SPRTN, match_any & has(RFLP_50)) # The RFLP number will change depending on the width
)
SPRTN
ggsave('figures/SPRTN.pdf', SPRTN, width = 8, height = 3)

# ---- SMARCAL1 ----
iSTOP_SMARCAL1 <-
  filter(CDS, gene == 'SMARCAL1') %>%
  locate_codons(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  locate_PAM(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  group_by(gene) %>%
  filter(cds_length == max(cds_length)) %>%
  add_RFLP(width = 50)

SMARCAL1 <- plot_spliced_isoforms(
  gene   = 'SMARCAL1',
  coords = filter(CDS, tx %in% iSTOP_SMARCAL1$tx),
  colors = c('black', 'blue', 'darkgreen'),
  #`Nonsense in cancer`   = COSMIC,
  `CAA, CAG, CGA, TGG`   = iSTOP_SMARCAL1,
  `iSTOP targetable`     = filter(iSTOP_SMARCAL1, match_any),
  `Verifiable with RFLP` = filter(iSTOP_SMARCAL1, match_any & has(RFLP_50)) # The RFLP number will change depending on the width
)
SMARCAL1
ggsave('figures/SMARCAL1.pdf', SMARCAL1, width = 8, height = 3)

# ---- CHEK2 ----
iSTOP_CHEK2 <-
  filter(CDS, gene == 'CHEK2') %>%
  locate_codons(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  locate_PAM(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  group_by(gene) %>%
  filter(cds_length == max(cds_length)) %>%
  add_RFLP(width = 50)

CHEK2 <- plot_spliced_isoforms(
  gene   = 'CHEK2',
  coords = filter(CDS, tx %in% iSTOP_CHEK2$tx),
  colors = c('black', 'blue', 'darkgreen'),
  #`Nonsense in cancer`   = COSMIC,
  `CAA, CAG, CGA, TGG`   = iSTOP_CHEK2,
  `iSTOP targetable`     = filter(iSTOP_CHEK2, match_any),
  `Verifiable with RFLP` = filter(iSTOP_CHEK2, match_any & has(RFLP_50)) # The RFLP number will change depending on the width
)
CHEK2
ggsave('figures/CHEK2.pdf', CHEK2, width = 8, height = 3)

# ---- PARP4 ----
iSTOP_PARP4 <-
  filter(CDS, gene == 'PARP4') %>%
  locate_codons(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  locate_PAM(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  group_by(gene) %>%
  filter(cds_length == max(cds_length)) %>%
  add_RFLP(width = 50)

PARP4 <- plot_spliced_isoforms(
  gene   = 'PARP4',
  coords = filter(CDS, tx %in% iSTOP_PARP4$tx),
  colors = c('black', 'blue', 'darkgreen'),
  #`Nonsense in cancer`   = COSMIC,
  `CAA, CAG, CGA, TGG`   = iSTOP_PARP4,
  `iSTOP targetable`     = filter(iSTOP_PARP4, match_any),
  `Verifiable with RFLP` = filter(iSTOP_PARP4, match_any & has(RFLP_50)) # The RFLP number will change depending on the width
)
PARP4
ggsave('figures/PARP4.pdf', PARP4, width = 8, height = 3)


# ---- PIK3R1 ----
iSTOP_PIK3R1 <-
  filter(CDS, gene == 'PIK3R1') %>%
  locate_codons(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  locate_PAM(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  group_by(gene) %>%
  filter(cds_length == max(cds_length)) %>%
  add_RFLP(width = 50)

PIK3R1 <- plot_spliced_isoforms(
  gene   = 'PIK3R1',
  coords = filter(CDS, tx %in% iSTOP_PIK3R1$tx),
  colors = c('black', 'blue', 'darkgreen'),
  #`Nonsense in cancer`   = COSMIC,
  `CAA, CAG, CGA, TGG`   = iSTOP_PIK3R1,
  `iSTOP targetable`     = filter(iSTOP_PIK3R1, match_any),
  `Verifiable with RFLP` = filter(iSTOP_PIK3R1, match_any & has(RFLP_50)) # The RFLP number will change depending on the width
)
PIK3R1
ggsave('figures/PIK3R1.pdf', PIK3R1, width = 8, height = 3)

# ---- ATM ----
iSTOP_ATM <-
  filter(CDS, gene == 'ATM') %>%
  locate_codons(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  locate_PAM(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  group_by(gene) %>%
  filter(cds_length == max(cds_length)) %>%
  add_RFLP(width = 50)

ATM <- plot_spliced_isoforms(
  gene   = 'ATM',
  coords = filter(CDS, tx %in% iSTOP_ATM$tx),
  colors = c('red', 'black', 'blue', 'darkgreen'),
  `Nonsense in cancer`   = COSMIC,
  `CAA, CAG, CGA, TGG`   = iSTOP_ATM,
  `iSTOP targetable`     = filter(iSTOP_ATM, match_any),
  `Verifiable with RFLP` = filter(iSTOP_ATM, match_any & has(RFLP_50)) # The RFLP number will change depending on the width
)
ATM
ggsave('figures/ATM.pdf', ATM, width = 8, height = 3)


# ---- SETD2 ----
iSTOP_SETD2 <-
  filter(CDS, gene == 'SETD2') %>%
  locate_codons(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  locate_PAM(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  group_by(gene) %>%
  filter(cds_length == max(cds_length)) %>%
  add_RFLP(width = 50)

SETD2 <- plot_spliced_isoforms(
  gene   = 'SETD2',
  coords = filter(CDS, tx %in% iSTOP_SETD2$tx),
  colors = c('red', 'black', 'blue', 'darkgreen'),
  `Nonsense in cancer`   = COSMIC,
  `CAA, CAG, CGA, TGG`   = iSTOP_SETD2,
  `iSTOP targetable`     = filter(iSTOP_SETD2, match_any),
  `Verifiable with RFLP` = filter(iSTOP_SETD2, match_any & has(RFLP_50)) # The RFLP number will change depending on the width
)
SETD2
ggsave('figures/SETD2.pdf', SETD2, width = 8, height = 3)

# ---- EZH2 ----
iSTOP_EZH2 <-
  filter(CDS, gene == 'EZH2') %>%
  locate_codons(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  locate_PAM(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  group_by(gene) %>%
  filter(cds_length == max(cds_length)) %>%
  add_RFLP(width = 50)

EZH2 <- plot_spliced_isoforms(
  gene   = 'EZH2',
  coords = filter(CDS, tx %in% iSTOP_EZH2$tx),
  colors = c('red', 'black', 'blue', 'darkgreen'),
  `Nonsense in cancer`   = COSMIC,
  `CAA, CAG, CGA, TGG`   = iSTOP_EZH2,
  `iSTOP targetable`     = filter(iSTOP_EZH2, match_any),
  `Verifiable with RFLP` = filter(iSTOP_EZH2, match_any & has(RFLP_50)) # The RFLP number will change depending on the width
)
EZH2
ggsave('figures/EZH2.pdf', EZH2, width = 8, height = 3)
