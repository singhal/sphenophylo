rm(list = ls())
library(sf)
library(ggplot2)
library(tidyverse)
library(ape)
library(phytools)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
source("scripts/colors.R")
cols = c("black", allcols[3], allcols[1])

#######

d = read.csv("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/metadata/sphenomorphine_all_v13.csv")
t = read.tree("data/diversification/species_level_phylogeny.OPERATIONAL_TAXON.tre")
x = read.csv("data/species_to_OTU/mapping.amended.csv")

# worst kluge to get missing support values
tx = read.tree("data/phylogenetic_support/synthetic_tree.all.support.dated.tre")
keep = d %>% filter(complete.cases(CURRENT_TAXON)) %>%
  group_by(CURRENT_TAXON) %>% slice_sample(n = 1) %>%
  pull(SAMPLE_ID)
tx = keep.tip(tx, keep)
tx$node.label = as.numeric(gsub("/\\S+", "", tx$node.label))
tx$tip.label = d[match(tx$tip.label, d$SAMPLE_ID), "OPERATIONAL_TAXON"]

for (i in 1:Nnode(tx)) {
  tips = tx$tip.label[ getDescendants(tx, Ntip(tx) + i)  ]
  tips1 = tips[complete.cases(tips)]
  
  if (length(tips1) < Ntip(t)) {
    t$node.label[ findMRCA(t, tips1, type = "node") - Ntip(t) ] = tx$node.label[i]
  }
}
t$node.label = as.numeric(t$node.label)


# get branch cols based on identity
bcols = rep("gray40", length(t$edge.length))
bcols[which(t$edge[,2] %in% match(x[which(x$delimitation == "same"), "OPERATIONAL_TAXON"], t$tip.label))] = cols[1]
bcols[which(t$edge[,2] %in% match(x[which(x$delimitation == "split"), "OPERATIONAL_TAXON"], t$tip.label))] = cols[2]
bcols[which(t$edge[,2] %in% match(x[which(x$delimitation == "combined"), "OPERATIONAL_TAXON"], t$tip.label))] = cols[3]

# make internal branches "lighter"
bwidth = rep(0.5, length(t$edge.length))
bwidth[which(t$edge[,2] <= Ntip(t))] = 0.8

pdf("figures/OTU_tree.pdf", height = 7, width = 6.5)
par(mar = c(2, 0, 0, 18), xpd = T)
plot.phylo(t, show.tip.label = F, edge.color = bcols, edge.width = bwidth)
lastPP = get("last_plot.phylo",envir=.PlotPhyloEnv)

# support values
nodelabels("", which(t$node.label >= 95) + Ntip(t), frame = "none", pch = 16, col = "black", cex = 0.3)

# genera names
for (gen in c("Ctenotus", "Lerista", "Notoscincus", "Concinnia", "Calyptotis", "Anomalopus",
              "Eulamprus", "Hemiergis", "Glaphyromorphus", "Eremiascincus")) {
  tips = grep(gen, t$tip.label)
  xloc1 = lastPP$x.lim[2] * 1.05
  lines(x = c(xloc1, xloc1), y = c(min(tips) + 0.5, max(tips) - 0.5))
  mid = mean(tips)
  xloc2 = lastPP$x.lim[2] * 1.07
  text(x = xloc2, y = mid, labels = gen, font = 3, adj = c(0, 0.5), cex = 0.7)
}

# legend
text(x = 1, y = 13, labels = "identical", 
     col = cols[1], adj = c(0, 0.5), cex = 0.7)
text(x = 1, y = 7, labels = "split", 
     col = cols[2], adj = c(0, 0.5), cex = 0.7)
text(x = 1, y = 1, labels = "combined", 
     col = cols[3], adj = c(0, 0.5), cex = 0.7)

# ages
maxht = max(phytools::nodeHeights(t))
# for (j in seq(5, 20, 5)) {
#  loc = max(lastPP$x.lim) - j
#  lines(x = c(loc, loc), y = c(-3, max(lastPP$y.lim)),
#        lty = "dotted")
# }
poly = c(seq(5, maxht, 5), round(maxht / 5, 0) * 5 + 5)
for (j in 1:length(poly)) {
  if (j %% 2 == 0) {
    loc1 = max(lastPP$x.lim) - (poly[j] - 5)
    loc2 = max(lastPP$x.lim) - poly[j]
    polygon(x = c(loc1, loc1, loc2, loc2), 
            y = c(-3, max(lastPP$y.lim), max(lastPP$y.lim), -3),
            col = alpha("black", 0.05),
            border = alpha("black", 0.05))
  }
}
 
xvals = c(0, 5, 10, 15, 20)
xvals2 = lastPP$x.lim[2] - xvals
axis(1, at = xvals2, labels = NA, line = 0.2, cex = 0.2, tck = -0.02)
axis(1, at = xvals2, labels = xvals, line= -0.9, lwd = 0, 
      cex.axis = 0.5)
mtext("time (Myr)", side=1, line=1, cex=0.7)
dev.off()
