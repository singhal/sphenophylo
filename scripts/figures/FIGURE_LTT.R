rm(list = ls())

library(ape)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(phytools)

get_ltt <- function(t) {
  t$node.label = seq(Ntip(t) + 1, Ntip(t) + Nnode(t))
  times = branching.times(t)
  
  tt = data.frame(
    node = t$node.label,
    age = times[as.character(t$node.label)]
  )
  tt[rev(order(tt$age)),  "ntip"] = seq(2, Ntip(t)) 
  maxht = max(nodeHeights(t))
  tt$time = -1 * tt$age
  tt = rbind(tt, c(NA, 0, Ntip(t), 0))
  return(tt)
}

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
source("scripts/colors.R")

d = read.csv("data/structure/inds_to_clusters.30May24.csv") %>% select(SAMPLE_ID = ind, pop)
x = read.csv("../Sphenomorphine_speciation/metadata/sphenomorphine_all_v11.csv")

t = read.tree("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")
outs = x[which(x$INGROUP == FALSE), "SAMPLE_ID"]
t = drop.tip(t, outs)

d1 = d[d$SAMPLE_ID %in% t$tip.label, ]
x1 = x[x$SAMPLE_ID %in% t$tip.label, ]

##############
# LTT all tips
##############

tt1 = get_ltt(t)
tt1$type = "clade-wide"

tt2 = get_ltt(keep.tip(t, t$tip.label[ which(x[match(t$tip.label, x$SAMPLE_ID), "GENUS"] %in% c("Lerista", "Ctenotus")) ]))
tt2$type = "Ctenotus & Lerista"

tt3 = get_ltt(keep.tip(t, t$tip.label[ which(!x[match(t$tip.label, x$SAMPLE_ID), "GENUS"] %in% c("Lerista", "Ctenotus")) ]))
tt3$type = "other genera"

tt = rbind(tt1, tt2, tt3)
a = ggplot(tt, aes(time, ntip)) + 
  geom_line(aes(col = type)) +
  theme_classic() +
  scale_y_log10() +
  xlab("time (myr before present)") +
  ylab("# of lineages") +
  scale_color_manual(values = allcols[1:3])

##############
# LTT - pops
##############

# first make subsampled tree
d1$species = gsub("_V\\d+", "", d1$pop)
inds1 = d1 %>% group_by(pop) %>%
  summarize(SAMPLE_ID = sample(SAMPLE_ID, 1)) %>% 
  pull(SAMPLE_ID)
# add in individuals not in clusters
inds2 = x1[!x1$CURRENT_TAXON %in% d1$species, ] %>% filter(complete.cases(CURRENT_TAXON)) %>% 
  group_by(CURRENT_TAXON) %>%
  summarize(SAMPLE_ID = sample(SAMPLE_ID, 1)) %>% 
  pull(SAMPLE_ID)
td = keep.tip(t, c(inds1, inds2))

bb1 = get_ltt(td)
bb1$type = "clade-wide"

bb2 = get_ltt(keep.tip(td, td$tip.label[ which(x[match(td$tip.label, x$SAMPLE_ID), "GENUS"] %in% c("Lerista", "Ctenotus")) ]))
bb2$type = "Ctenotus & Lerista"

bb3 = get_ltt(keep.tip(td, td$tip.label[ which(!x[match(td$tip.label, x$SAMPLE_ID), "GENUS"] %in% c("Lerista", "Ctenotus")) ]))
bb3$type = "other genera"

bb = rbind(bb1, bb2, bb3)
b = ggplot(bb, aes(time, ntip)) + 
  geom_line(aes(col = type)) +
  theme_classic() +
  scale_y_log10() +
  xlab("time (myr before present)") +
  ylab("# of lineages") +
  scale_color_manual(values = allcols[1:3])

##############
# LTT - OTUs
##############

tc = read.tree("data/diversification/species_level_phylogeny.CURRENT_TAXON.tre")

aa1 = get_ltt(tc)
aa1$type = "clade-wide"

aa2 = get_ltt(keep.tip(tc, tc$tip.label[c(grep("Ctenotus", tc$tip.label), grep("Lerista", tc$tip.label))]))
aa2$type = "Ctenotus & Lerista"

aa3 = get_ltt(keep.tip(tc, tc$tip.label[-c(grep("Ctenotus", tc$tip.label), grep("Lerista", tc$tip.label))]))
aa3$type = "other genera"

aa = rbind(aa1, aa2, aa3)
c = ggplot(aa, aes(time, ntip)) + 
  geom_line(aes(col = type))  +
  scale_y_log10() +
  xlab("time (myr before present)") +
  ylab("# of lineages") +
  scale_color_manual(values = allcols[1:3]) +
  theme_classic() +
  theme(legend.position = c(0.3, 0.87),
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA))
  
prow <- plot_grid(
  a + theme(legend.position="none"),
  b + theme(legend.position="none"),
  c,
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)
save_plot("figures/LTT.png", prow, base_height = 3, base_width = 10)
