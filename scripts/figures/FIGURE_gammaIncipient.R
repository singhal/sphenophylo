rm(list = ls())

library(ape)
library(ggplot2)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")

source("scripts/colors.R")

res2 = read.csv("data/gamma/gamma_incipient_sim.csv")

phy = read.tree("data/diversification/species_level_phylogeny.CLUSTERS.tre")
gamma = gammaStat(phy)

phy = read.tree("data/diversification/species_level_phylogeny.OPERATIONAL_TAXON.tre")
gamma2 = gammaStat(phy)

ab = ggplot(res2, aes(as.factor(X1), X3)) + 
  geom_boxplot() +
  geom_hline(aes(yintercept = gamma), linetype = 2, color = tcols["incipient"]) +
  annotate("text", y=gamma + 0.7, x=25, hjust = 1, label= "incipient", col = tcols["incipient"]) + 
  geom_hline(aes(yintercept = gamma2), linetype = 2, color = tcols["operational"]) +
  annotate("text", y=gamma2 + 0.7, x=25, hjust = 1, label= "operational", col = tcols["operational"]) +
  xlab("# of incipient lineages") +
  ylab(expression(gamma)) +
  theme_bw()
cowplot::save_plot("figures/gamma_incipient.png", ab, 
                   base_height = 3, base_width = 7)
