rm(list = ls())
library(sf)
library(ggplot2)
library(tidyverse)
library(ape)
library(phytools)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
source("scripts/colors.R")

#######

d = read.csv("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/metadata/sphenomorphine_all_v10.csv")
t = read.tree("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")

outs = d[which(d$INGROUP == FALSE), "SAMPLE_ID"]
t = drop.tip(t, outs)
t$node.label = seq(Ntip(t) + 1, Ntip(t) + Nnode(t))
times = branching.times(t)

tt = data.frame(
  node = t$node.label,
  age = times[as.character(t$node.label)]
)
tt[rev(order(tt$age)),  "ntip"] = seq(2, Ntip(t)) 
maxht = max(nodeHeights(t))
tt$time = -1 * tt$age

tt$type = "all"
# identify within
for (i in 1:nrow(tt)) {
  node = tt$node[i]
  desc = phytools::getDescendants(t, node)
  desc = t$tip.label[ desc[desc < Ntip(t) + 1] ]
  sps = unique(d[d$SAMPLE_ID %in% desc, "CURRENT_TAXON"])
  if (length(sps) == 1) {
    tt[tt$node == node, "type"] = "within"
  }
}
sps = unique(d$CURRENT_TAXON)
sps = sps[complete.cases(sps)]
# identify crown ages & stem (sister) ages
for (sp in sps) {
  desc = d[which(d$CURRENT_TAXON == sp), "SAMPLE_ID"]
  if (length(desc) == 1) {
    parent = phytools::getParent(t, which(t$tip.label == desc))
    tt[tt$node == parent, "type"] = "stem"
  } else {
    node = phytools::findMRCA(t, desc)
    parent = phytools::getParent(t, node)
    # cat(sp, "\n")
    tt[tt$node == node, "type"] = "crown"
    tt[tt$node == parent, "type"] = "stem"
  }
}

table(tt$type)

a = ggplot(tt, aes(time, ntip)) +
  geom_point(col = "gray80")  +
  geom_point(data = tt %>% filter(type == "within"), 
             aes(time, ntip),
             fill = genera[3], pch = 21)  +
  geom_point(data = tt %>% filter(type == "within"), 
             aes(time, ntip),
             col = genera[3], pch = 16)  +
  xlab("time (myr)") +
  ylab("# of tips") +
  theme_classic()

names(genera) = NULL
b = ggplot(tt %>% filter(type != "all"), aes(age)) + 
  geom_density(aes(fill = type), alpha = 0.3, linewidth = 0.1) +
  xlab("age (Myr)") + 
  scale_fill_manual(values = genera) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.75))

ab = cowplot::plot_grid(a, b, labels = c("A", "B"))
cowplot::save_plot("figures/breaks.png", ab, base_width = 6.5, base_height = 2.5)
