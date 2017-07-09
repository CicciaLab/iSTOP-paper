species_order <- c(
  'Hsapiens'      = 'H. sapiens',
  'Mmusculus'     = 'M. musculus',
  'Rnorvegicus'   = 'R. norvegicus',
  'Drerio'        = 'D. rerio',
  'Dmelanogaster' = 'D. melanogaster',
  'Celegans'      = 'C. elegans',
  'Athaliana'     = 'A. thaliana',
  'Scerevisiae'   = 'S. cerevisiae'
)

codon <- read_csv('data/Figure-data/Codon-Summary.csv') %>% mutate(species = ordered(species, levels = names(species_order), labels = species_order))
tx    <- read_csv('data/Figure-data/ORF-Summary.csv') %>% mutate(species = ordered(species, levels = names(species_order), labels = species_order))
gene  <- read_csv('data/Figure-data/Gene-Summary.csv') %>%mutate(species = ordered(species, levels = names(species_order), labels = species_order))

# ---- Codons ----
codon_summary_all_species <-
  codon %>%
  select(species, Available = Targetable, Unavailable = `Not Targetable`) %>%
  mutate(species = forcats::fct_rev(species)) %>%
  gather(`PAM:`, n, -species) %>%
  mutate(`PAM:` = factor(`PAM:`, levels = c('Unavailable', 'Available'))) %>%
  ggplot(aes(y = n, x = species, fill = `PAM:`)) +
  geom_col(color = 'black', position = 'stack', width = 0.6) +
  geom_text(aes(x = species, y = Total,        label = paste0(round((`Not Targetable` / Total) * 100, digits = 1), '%')), data = codon, inherit.aes = F, color = 'black', fontface = 'bold', hjust = 1.1) +
  geom_text(aes(x = species, y = `Targetable`, label = paste0(round((Targetable / Total) * 100,       digits = 1), '%')), data = codon, inherit.aes = F, color = 'white', fontface = 'bold', hjust = 1.1) +
  scale_fill_manual(values = c('white', '#77cd99'), labels = c('Untargetable', 'Targetable')) +
  scale_y_continuous(breaks = seq(0, 800000, by = 100000), expand = c(0, 0), labels = scales::comma) +
  coord_flip(ylim = c(0, 900000)) +
  labs(x = NULL, y = 'CAA, CAG, CGA, and TGG codons') +
  theme_bw(base_size = 16) +
  guides(fill = guide_legend(title = NULL, reverse = T)) +
  theme(
    axis.text.y = element_text(face = 'italic', color = 'black', hjust = 0.5, vjust = 0.5),
    plot.margin = unit(c(1,0.25,0.25,0),"cm"),
    panel.grid.major.x = element_line(color = 'black', linetype = 'dotted'),
    legend.title = element_text(face = 'bold'),
    legend.direction = 'vertical',
    legend.position = c(1, 0),
    legend.justification = c(1.1, -0.1),
    legend.background = element_rect(colour = 'black', fill = 'white', size = 0.25)
  )

codon_summary_all_species

#gg_aa_per_site <-
#  codon %>%
#  select(species, sites_per_100aa) %>%
#  mutate(
#    species = forcats::fct_rev(species),
#    aa_per_site = 1 / (sites_per_100aa / 100)
#  ) %>%
#  ggplot(aes(y = aa_per_site, x = species)) +
#  geom_col(color = 'black', fill = '#77cd99', width = 0.6) +
#  geom_text(aes(x = species, y = aa_per_site, label = round(aa_per_site, digits = 1)), inherit.aes = F, color = 'white', fontface = 'bold', hjust = -0.2) +
#  scale_y_reverse(breaks = seq(0, 35, by = 10), expand = c(0, 0)) +
#  scale_x_discrete(position = 'top') +
#  coord_flip(ylim = c(0, 35)) +
#  theme_bw(base_size = 16) +
#  labs(x = NULL, y = 'Codons per targetable site') +
#  theme(
#    #axis.text.y = element_text(face = 'italic', color = 'black', hjust = 0.5, vjust = 0.5),
#    axis.text.y = element_blank(),
#    plot.margin = unit(c(1,0.1,0.25,0.25),"cm"),
#    #panel.grid.minor = element_blank(),
#    panel.grid.major.x = element_line(color = 'black', linetype = 'dotted')
#  )
#gg_aa_per_site

#gg1 <-
#  ggdraw() +
#  draw_plot(gg_aa_per_site, x = 0,    y = 0, width = 0.25, height = 1) +
#  draw_plot(gg_codons,         x = 0.25, y = 0, width = 0.75, height = 1) +
#  draw_plot_label(c('A')) +
#  draw_plot_label('Codons', x = 0.03, fontface = 'bold', hjust = -0.1, vjust = 2, size = 12)
#gg1

ggsave('figures/Codon-summary-all-species.pdf', codon_summary_all_species, height = 4, width = 11)


# ---- ORFs ----
ORF_summary_all_species_targetable <-
  tx %>%
  select(species, Targetable, `% Targetable`) %>%
  mutate(species = forcats::fct_rev(species)) %>%
  ggplot(aes(y = Targetable, x = species)) +
  geom_col(color = 'black', fill = '#8aa6da', position = 'stack', width = 0.6) +
  geom_text(aes(label = paste0(round(`% Targetable`, digits = 1), '%')), color = 'white', fontface = 'bold', hjust = 1.1) +
  scale_y_continuous(breaks = seq(0, 70000, by = 10000), expand = c(0, 0), labels = scales::comma) +
  coord_flip(ylim = c(0, 75000)) +
  labs(x = NULL, y = 'Targetable ORFs') +
  theme_bw(base_size = 16) +
  theme(
    axis.text.y = element_text(face = 'italic', color = 'black', hjust = 0.5, vjust = 0.5),
    plot.margin = unit(c(1,0.25,0.25,0),"cm"),
    #panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = 'black', linetype = 'dotted'),
    legend.title = element_blank(),
    legend.position = c(1, 0),
    legend.justification = c(1.1, -0.1),
    legend.background = element_rect(colour = 'black', fill = 'white', size = 0.25)
    #    axis.ticks.length = unit(0.3, 'cm')
  )
ORF_summary_all_species_targetable


ORF_summary_all_species_untargetable <-
  tx %>%
  select(species, `PAM` = `Not Targetable`, `Codon` = `No Codons`) %>%
  gather(type, n, -species) %>%
  mutate(species = forcats::fct_rev(species)) %>%
  ggplot(aes(y = n, x = species, fill = type)) +
  geom_col(color = 'black', position = 'identity', width = 0.6) +
  geom_text(aes(x = species, y = `Not Targetable`, label = paste0(round(((`Not Targetable` - `No Codons`) / `Total`) * 100, digits = 1), '%')), inherit.aes = F, color = 'black', fontface = 'bold', hjust = -0.2, data = tx) +
  #geom_text(aes(x = species, y = `No Codons`, label = paste0(round((`No Codons` / `Total`) * 100, digits = 1), '%')), inherit.aes = F, color = 'white', fontface = 'bold', hjust = -0.2, data = tx) +
  scale_y_continuous(breaks = seq(0, 1000, by = 200), expand = c(0, 0), labels = scales::comma) +
  #scale_x_discrete(position = 'top') +
  scale_fill_manual(values = c('grey', 'white')) +
  coord_flip(ylim = c(0, 1100)) +
  labs(x = NULL, y = 'Untargetable') +
  guides(fill = guide_legend(title = 'Due to unavailable:')) +
  theme_bw(base_size = 16) +
  theme(
    #axis.text.y = element_blank(),
    #axis.ticks.y = element_blank(),
    plot.margin = unit(c(1,0.1,0.25,0.25),"cm"),
    #panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = 'black', linetype = 'dotted'),
    legend.position = c(0.5, 0.5),
    legend.justification = c(-0.1, -0.1),
    legend.background = element_rect(colour = 'black', fill = 'white', size = 0.25)
    #    axis.ticks.length = unit(0.3, 'cm')
  )
ORF_summary_all_species_untargetable


#gg2 <-
#  ggdraw() +
#  draw_plot(gg_transcripts_untargetable, x = 0,    y = 0, width = 0.25, height = 1) +
#  draw_plot(gg_transcripts_targetable,   x = 0.25, y = 0, width = 0.75, height = 1) +
#  draw_plot_label(c('B')) +
#  draw_plot_label('ORFs', x = 0.03, fontface = 'bold', hjust = -0.1, vjust = 2, size = 12)
#gg2

ggsave('figures/ORF-summary-all-species-targetable.pdf', ORF_summary_all_species_targetable, height = 4, width = 11)
ggsave('figures/ORF-summary-all-species-untargetable.pdf', ORF_summary_all_species_untargetable, height = 5, width = 5)


# ---- Genes ----
genes_targetable_all_species <-
  gene %>%
  select(species, one_or_more_targetable, all_targetable) %>%
  gather(type, n, -species) %>%
  mutate(species = forcats::fct_rev(species), type = if_else(type == 'all_targetable', 'All', '1+')) %>%
  ggplot(aes(y = n, x = species, fill = type)) +
  geom_col(color = 'black', position = 'identity', width = 0.6) +
  geom_col(aes(y=0, fill="0"),size=3) +
  geom_text(aes(x = species, y = all_targetable, label = paste0(round((all_targetable / Total) * 100, digits = 1), '%')), inherit.aes = F, color = 'black', fontface = 'bold', hjust = 1.2, data = gene) +
  geom_text(aes(x = species, y = one_or_more_targetable, label = paste0(round(((one_or_more_targetable - all_targetable) / Total) * 100, digits = 1), '%')), inherit.aes = F, color = 'black', fontface = 'bold', hjust = -0.2, data = gene) +
  scale_y_continuous(breaks = seq(0, 20000, by = 2500), expand = c(0, 0), labels = scales::comma) +
  scale_fill_manual(values = c('grey', '#d02c1e', '#ee8671')) +
  coord_flip(ylim = c(0, 22000)) +
  labs(x = NULL, y = 'Targetable genes') +
  guides(fill = guide_legend(title = 'Isoforms:', reverse = T)) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.y = element_text(face = 'italic', color = 'black', hjust = 0.5, vjust = 0.5),
    plot.margin = unit(c(1,0.25,0.25,0),"cm"),
    #panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = 'black', linetype = 'dotted'),
    legend.position = c(1, 0),
    legend.justification = c(1.1, -0.1),
    legend.direction = 'horizontal',
    legend.background = element_rect(colour = 'black', fill = 'white', size = 0.25)
    #    axis.ticks.length = unit(0.3, 'cm')
  )
genes_targetable_all_species


#gg_genes_untargetable <-
#  gene %>%
#  select(species, Total, one_or_more_targetable) %>%
#  mutate(species = forcats::fct_rev(species), `Not Targetable` = Total - one_or_more_targetable) %>%
#  ggplot(aes(y = `Not Targetable`, x = species)) +
#  geom_col(color = 'black', fill = 'grey', position = 'identity', width = 0.6) +
#  geom_text(aes(label = paste0(round((`Not Targetable` / `Total`) * 100, digits = 1), '%')), color = 'black', fontface = 'bold', hjust = -0.2) +
#  scale_y_reverse(breaks = seq(0, 200, by = 50), expand = c(0, 0), labels = scales::comma) +
#  scale_x_discrete(position = 'top') +
#  #scale_fill_manual(values = c('grey', 'white')) +
#  coord_flip(ylim = c(0, 250)) +
#  labs(x = NULL, y = 'Untargetable') +
#  guides(fill = guide_legend(title = 'Due to unavailable:')) +
#  theme_bw(base_size = 16) +
#  theme(
#    axis.text.y = element_blank(),
#    plot.margin = unit(c(1,0.1,0.25,0.25),"cm"),
#    #panel.grid.minor = element_blank(),
#    panel.grid.major.x = element_line(color = 'black', linetype = 'dotted'),
#    legend.position = c(1, 0.5),
#    legend.justification = c(1.1, -0.1),
#    legend.background = element_rect(colour = 'black', fill = 'white', size = 0.25)
#    #    axis.ticks.length = unit(0.3, 'cm')
#  )
#gg_genes_untargetable

#gg3 <-
#  ggdraw() +
#  draw_plot(gg_genes_untargetable, x = 0,    y = 0, width = 0.25, height = 1) +
#  draw_plot(gg_genes_targetable,   x = 0.25, y = 0, width = 0.75, height = 1) +
#  draw_plot_label(c('C')) +
#  draw_plot_label('Genes', x = 0.03, fontface = 'bold', hjust = -0.1, vjust = 2, size = 12)
#gg3

ggsave('figures/Genes-summary-all-species-targetable.pdf', genes_targetable_all_species, height = 4, width = 11)

# ---- Combined ----

#combined <-
#  ggdraw() +
#  draw_plot(gg1, x = 0, y = 0.666, width = 1, height = 0.333) +
#  draw_plot(gg2, x = 0, y = 0.333, width = 1, height = 0.333) +
#  draw_plot(gg3, x = 0, y = 0,     width = 1, height = 0.333)
#combined
#
#ggsave('figures/Figure-4.pdf', combined, height = 11 * 1.5, width = 8.5 * 1.5)

