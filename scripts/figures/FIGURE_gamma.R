rm(list = ls())

library(tidyverse)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")

source("scripts/colors.R")

res2 = readRDS("~/Desktop/gamma_simulations.Rds")
res3 = res2 %>% group_by(agedrop_per, agedrop) %>% 
  summarize(n = n()) %>% 
  mutate(drop = ifelse(n < 20, TRUE, FALSE))
res4 = res2 %>% left_join(res3) %>% filter(drop == FALSE)

t2 = read.tree("data/diversification/species_level_phylogeny.OPERATIONAL_TAXON.tre")
t1 = read.tree("data/diversification/species_level_phylogeny.MORPHOLOGICAL_TAXON.tre")
t3 = read.tree("data/diversification/species_level_phylogeny.CLUSTERS.tre")
t4 = read.tree("data/diversification/species_level_phylogeny.THRESHOLD_TAXON.tre")

vals = c("strong", "weak", "random")
names(vals) = c(0.25, 0.5, 1)

res4$agedrop2 = factor(vals[ as.character(res4$agedrop) ], levels = vals)

x = read.csv("data/operational_taxonomy/operational_taxonomy.csv")
persample = Ntip(t2) / nrow(x)

a = ggplot(res4, aes(as.factor(agedrop_per * 100), agedrop_tree_gamma, 
                 fill = agedrop2)) +
  geom_boxplot(outliers = FALSE, lwd = 0.3) + 
  ylab(expression(gamma~"-statistic")) +
  geom_hline(aes(yintercept = gammaStat(t1)), linetype = 2, color = tcols["morphological"]) + 
  annotate("text", y=gammaStat(t1) + 0.7, x=11, hjust = 1, label= "morphological", col = tcols["morphological"]) + 
  geom_hline(aes(yintercept = gammaStat(t2)), linetype = 2, color = tcols["operational"]) + 
  annotate("text", y=gammaStat(t2) + 0.7, x=11, hjust = 1,  label= "operational", col = tcols["operational"]) + 
  geom_hline(aes(yintercept = gammaStat(t3)), linetype = 2, color = tcols["incipient"]) + 
  annotate("text", y=gammaStat(t3) + 0.7, x=11, hjust = 1,  label= "incipient", col = tcols["incipient"]) +
  geom_hline(aes(yintercept = gammaStat(t4)), linetype = 2, color = tcols["threshold"]) + 
  annotate("text", y=gammaStat(t4) - 0.7, x=11, hjust = 1,  label= "threshold", col = tcols["threshold"]) +
  xlab("% of species sampled") + 
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Greys")) +
  scale_colour_manual(values = RColorBrewer::brewer.pal(5, "Greys")) +
  theme_classic() +
  guides(fill=guide_legend(title="bias in dropped\nlineages"))

cowplot::save_plot("figures/gamma.png", a, base_width = 9, base_height = 3)
