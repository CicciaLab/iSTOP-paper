#library(igraph)
#library(ggraph)

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
ggsave('figures/Supplemental-Figure-TSG-Heatmap.pdf', Supp_fig_TSG_Heatmap, width = 16, height = 3)

with(Comparison, table(TSG_iSTOP,  TSG_TUSCON))
with(Comparison, table(TSG_iSTOP,  TSG_20_20))
with(Comparison, table(TSG_iSTOP,  TSG_CGC))
with(Comparison, table(TSG_TUSCON, TSG_CGC))
with(Comparison, table(TSG_20_20,  TSG_CGC))
with(Comparison, table(TSG_TUSCON, TSG_20_20))



iSTOP_frequent <-
  COSMIC_by_gene %>%
  filter(q < 0.05 | cancer_type == 'All cancers') %>%
  group_by(gene) %>%
  summarise(
    n_cancers                                = sum(q[which(cancer_type != 'All cancers')] < 0.05),
    q                                        = q[which(cancer_type == 'All cancers')], # strongest score
    `-log10(q)`                              = -log10(q),
    cancer_types                             = str_c(c(cancer_type[which(q < 0.05 & cancer_type != 'All cancers')], ''), collapse = ' | '),
    n_nonsense_in_gene_in_cancer             = n_nonsense_in_gene_in_cancer[cancer_type == 'All cancers'],
    n_iSTOP_sites_in_gene                    = n_iSTOP_sites_in_gene[cancer_type == 'All cancers'],
    n_iSTOP_sites_in_gene_in_cancer          = n_iSTOP_sites_in_gene_in_cancer[cancer_type == 'All cancers'],
    n_iSTOP_sites_in_gene_in_cancer_with_PAM = n_iSTOP_sites_in_gene_in_cancer_with_PAM[cancer_type == 'All cancers'],
    frac_iSTOP_in_cancer                     = n_iSTOP_sites_in_gene_in_cancer / n_iSTOP_sites_in_gene,
    frac_iSTOP_in_cancer_with_PAM            = n_iSTOP_sites_in_gene_in_cancer_with_PAM / n_iSTOP_sites_in_gene
  ) %>%
  mutate(cancer_types = str_replace(cancer_types, ' [|] $', '')) %>%
  filter(!cancer_types %in% c('All cancers', ''))


iSTOP_frequent_gene <-
  COSMIC_by_gene %>%
  filter(q < 0.05 | cancer_type == 'All cancers') %>%
  group_by(gene) %>%
  mutate(n_cancers = sum(q < 0.05)) %>%
  filter(n_cancers > 0, !(n_cancers == 1 & cancer_type == 'All cancers'))


simple_roc <- function(labels, scores){
  scores_order <- order(scores, decreasing=TRUE)
  labels <- labels[scores_order]
  scores <- scores[scores_order]
  data.frame(
    TPR = cumsum(labels)  / sum(labels),
    FPR = cumsum(!labels) / sum(!labels),
    labels,
    scores
  )
}

x <- filter(iSTOP_frequent, `-log10(q)` > 3) %>% mutate(TSG = gene %in% filter(CGC, str_detect(`Role in Cancer`, 'TSG'))$`Gene Symbol`) %>% select(TSG, everything())
nrow(x)
nrow(filter(CGC, str_detect(`Role in Cancer`, 'TSG')))
(nrow(x) - sum(!x$TSG)) / nrow(x)   # 1/3 are classified as TSG. This captures 40 TSG of a possible 132 (30% of TSG)


df <-
  COSMIC_by_gene %>%
  filter(cancer_type != 'All cancers') %>%
  group_by(gene) %>%
  summarise(score = max(`-log10(q)`)) %>%
  mutate(TSG = gene %in% filter(CGC, str_detect(`Role in Cancer`, 'TSG'))$`Gene Symbol`) %>%
  ungroup %>%
  arrange(desc(score))

simple_roc(df$TSG, df$score) %>%
  ggplot(aes(x = FPR, y = TPR, color = scores > 3)) +
  geom_line(aes(alpha = scores > 3, linetype = scores > 3), size = 1) +
  geom_segment(aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = F, data = data_frame(x = 0, y = 0, xend = 1, yend = 1)) +
  labs(x = 'False Positive Rate', y = 'True Positive Rate', title = 'ROC curve', subtitle = 'True positives: Cancer Gene Census tumor suppressor genes (TSG)') +
  theme_bw() +
  scale_linetype_manual(values = c('dotted', 'solid')) +
  scale_alpha_discrete(range = c(0.6, 1)) +
  guides(color = guide_legend(title = 'Frequent iSTOP\n-log10(q) > 3'), alpha = guide_legend(title = 'Frequent iSTOP\n-log10(q) > 3'), linetype = guide_legend(title = 'Frequent iSTOP\n-log10(q) > 3')) +
  coord_fixed() +
  theme(legend.position = c(0.95, 0.05), legend.justification = c(1, 0), legend.background = element_rect(color = 'black'))

ggsave('figures/ROC-frequent-iSTOPers.pdf', width = 5, height = 5)
