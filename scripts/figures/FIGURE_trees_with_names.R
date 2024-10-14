rm(list = ls())
library(ape)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
source("scripts/colors.R")

#######

d = read.csv("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/metadata/sphenomorphine_all_v12.csv")
ct = read.tree("data/diversification/species_level_phylogeny.OPERATIONAL_TAXON.tre")
nt = read.tree("data/diversification/species_level_phylogeny.MORPHOLOGICAL_TAXON.tre")

pdf("figures/trees_with_name.pdf", height = 12, width = 6.5)
par(mfrow = c(1, 2), mar = c(0, 0, 0, 8), xpd = T)
plot.phylo(ct, show.tip.label = F, edge.width = 0.5)
tiplabels(gsub("_", " ", ct$tip.label), frame = "none", adj = c(0, 0.5), cex = 0.4, font = 3)
plot.phylo(nt, show.tip.label = F, edge.width = 0.5)
tiplabels(gsub("_", " ", nt$tip.label), frame = "none", adj = c(0, 0.5), cex = 0.5, font = 3)
dev.off()