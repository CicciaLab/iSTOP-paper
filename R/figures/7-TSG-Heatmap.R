Twenty_Twenty <-
  readxl::read_excel('data/PMID27911828/pnas.1616440113.sd04.xlsx', sheet = 1, skip = 1) %>%
  select(gene, q_TSG_20_20 = `tsg q-value`)

TUSON <-
  readxl::read_excel('data/PMID27911828/pnas.1616440113.sd04.xlsx', sheet = 2, skip = 1) %>%
  select(gene = Gene, q_TSG_TUSON = TUSON.combined.qvalue.TSG)

MutSigCV <-
  readxl::read_excel('data/PMID27911828/pnas.1616440113.sd04.xlsx', sheet = 3, skip = 1) %>%
  select(gene, q_MutSig = q)

TSG <- list(Twenty_Twenty, TUSON, MutSigCV) %>% reduce(inner_join) #inner_join(Twenty_Twenty, TUSON)
CGC <- read_csv('data/COSMIC/CGC.csv')
COSMIC_by_gene <- read_csv('data/COSMIC/COSMIC-summary-by-gene.csv')

Comparison <-
  COSMIC_by_gene %>%
  filter(cancer_type == 'All cancers') %>%
  select(gene, q) %>%
  inner_join(TSG) %>%
  mutate(
    score_iSTOP    = -log10(q),
    score_TUSCON   = -log10(q_TSG_TUSON),
    score_20_20    = -log10(q_TSG_20_20),
    score_MutSigCV = -log10(q_MutSig),
    TSG_iSTOP      = score_iSTOP    >= sort(score_iSTOP,  decreasing = T)[100],
    TSG_TUSON      = score_TUSCON   >= sort(score_TUSCON, decreasing = T)[100],
    TSG_20_20      = score_20_20    >= sort(score_20_20,  decreasing = T)[100],
    MutSigCV       = score_MutSigCV >= sort(score_MutSigCV, decreasing = T)[100],
    TSG_CGC        = gene %in% filter(CGC, str_detect(`Role in Cancer`, 'TSG'))$`Gene Symbol`
  )


Supp_Table_TSG <-
  Comparison %>%
  filter(TSG_iSTOP | TSG_TUSON | TSG_20_20 | MutSigCV) %>%
  write_csv('data/Figure-data/Supp-Table-TSG.csv')


TSG_mat <-
  Supp_Table_TSG %>%
  mutate(n = TSG_CGC + TSG_iSTOP + TSG_20_20 + TSG_TUSON + MutSigCV) %>%
  arrange(-TSG_CGC, -n, -TSG_iSTOP, -TSG_20_20, -TSG_TUSON, -MutSigCV, score_iSTOP) %>%
  select(gene, starts_with('TSG'), MutSigCV, n) %>%
  filter(n > 2 | TSG_CGC)

Supp_fig_TSG_Heatmap <-
  Supp_Table_TSG %>%
  select(gene, starts_with('TSG'), MutSigCV) %>%
  gather(key, value, -gene) %>%
  filter(gene %in% TSG_mat$gene) %>%
  mutate(
    key  = ordered(key,  levels = c('TSG_CGC', 'TSG_iSTOP', 'TSG_20_20', 'TSG_TUSON', 'MutSigCV', 'n'), labels = c('Cancer Gene Census', 'iSTOPers', '20/20+', 'TUSON', 'MutSigCV', 'Total')),
    gene = ordered(gene, levels = TSG_mat$gene),
    key  = forcats::fct_rev(key)
  ) %>%
  ggplot(aes(x = gene, y = key, fill = value)) +
  geom_tile(color = 'white') +
  scale_fill_manual(values = c('white', 'red')) +
  #coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank())
Supp_fig_TSG_Heatmap
ggsave('figures/TSG-Heatmap.pdf', Supp_fig_TSG_Heatmap, width = 16, height = 3)

