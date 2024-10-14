rm(list = ls())

library(ape)
library(diversitree)
library(tidyverse)
library(phytools)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")

calculate_gamma_sim <- function(ntip, miss, agedrop) {
  # simulate the tree
  total_tip = round(ntip / miss)
  tree <- tree.bd(pars = c(1, 0), max.taxa = total_tip,
                  max.t = Inf, include.extinct = FALSE)
  
  # gamma on the completely sampled tree
  gamma1 <- gammaStat(tree)
  
  # gamma on the partially sampled tree
  dropn = total_tip - ntip
  tree2 = drop.tip(tree, tip = sample(tree$tip.label, dropn))
  gamma2 = gammaStat(tree2)
  
  # gamma on the tree with just recent tips unsampled
  hts = nodeHeights(tree)
  tipname = tree$tip.label[ tree$edge[,2][ which(tree$edge[,2] <= Ntip(tree)) ] ]
  tiphts = data.frame(hts[which(tree$edge[,2] <= Ntip(tree)), ],
                      tipname)
  tiphts$age = tiphts$X2 - tiphts$X1
  young = tiphts[tiphts$age <= agedrop * max(tiphts$age), ]
  if (dropn > nrow(young)) {
    dropn = nrow(young)
  }
  tree3 = drop.tip(tree, tip = sample(young$tipname, dropn))
  gamma3 = gammaStat(tree3)
  
  res = c(Ntip(tree), gamma1, Ntip(tree2), miss, gamma2, Ntip(tree3), agedrop, gamma3)
  names(res) = c("full_tree_n", "full_tree_gamma",
                 "subsample_tree_n", "subsample_per", "subsample_tree_gamma",
                 "agedrop_tree_n", "agedrop", "agedrop_tree_gamma")
  return(res)
}

# starting trees
t1 = read.tree("data/diversification/species_level_phylogeny.MORPHOLOGICAL_TAXON.tre")
t2 = read.tree("data/diversification/species_level_phylogeny.OPERATIONAL_TAXON.tre")
t3 = read.tree("data/diversification/species_level_phylogeny.CLUSTERS.tre")

# percentage of tips sampled
missc = seq(0.1, 1, 0.1)

# num of simulations
nsim = 1000

# agedrop - relative to oldest split, not total tree height
agedropc = c(0.25, 0.5, 1)

# store results
res = vector("list", nsim * length(missc) * length(agedropc))
x = 1
for (miss in missc) {
  for (agedrop in agedropc) {
    message(paste0("doing ", miss, " & ", agedrop))
    for (i in 1:nsim) {
      res[[x]] = calculate_gamma_sim(Ntip(t2), miss, agedrop)
      x = x + 1
    }
  }
}

res2 = data.frame(do.call("rbind", res))
res2$agedrop_per = round(res2$agedrop_tree_n / res2$full_tree_n, 1)

saveRDS(res2, "~/Desktop/gamma_simulations.Rds")
