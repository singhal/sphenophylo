rm(list = ls())

library(ape)
library(diversitree)
library(tidyverse)
library(phytools)
library(cowplot)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")

source("scripts/DR_functions.R")

get_groups <- function(t, breakpt) {
  # http://blog.phytools.org/2014/11/pruning-trees-to-one-member-per-genus.html
  H = nodeHeights(t)
  threshold = max(H) - breakpt
  h1<-which(H[,1]<threshold)
  h2<-which(H[,2]>threshold)
  ii<-intersect(h1,h2)
  ## all daughter nodes of those edges
  nodes<-t$edge[ii,2]
  getDescendants<-phytools:::getDescendants
  ## find all descendants from each edge
  tips <- lapply(nodes,getDescendants,tree=t)
  tips2 = lapply(tips, function(x) as.character(na.omit(t$tip.label[x])))
  names(tips2) = paste0("sp", seq(1:length(tips2)))
  tips3 = data.frame(sp = unlist(mapply(function(x, y) {rep(x, length(y))}, names(tips2), tips2)),
                     ind = unlist(tips2))
  return(tips3)
}

t = read.tree("../Sphenomorphine_speciation/phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")
d = read.csv("../Sphenomorphine_speciation/metadata/sphenomorphine_all_v13.csv")

outs = d[which(d$INGROUP == FALSE), "SAMPLE_ID"]
t = drop.tip(t, outs)

thresh = seq(0.1, 10, 0.1)
res = vector("list", length(thresh))
res2 = vector("list", length(thresh))
# threshold the tree at different levels
# calculate DR and gamma
for (i in 1:length(thresh)) {
  gps = get_groups(t, thresh[i])
  # subsample tree
  inds = gps %>% group_by(sp) %>% 
    summarize(sample = sample(ind, 1)) %>% pull(sample)
  t1 = keep.tip(t, inds)
  
  # gamma stat
  gamma = gammaStat(t1)
  
  # DR statistic 
  r1 = data.frame(DR(t1)) 
  
  rr = right_join(gps, r1 %>% mutate(ind = row.names(r1)))
  rr1 = rr$DR.t1.
  names(rr1) = rr$sp
  gps$rates = rr1[ gps$sp ]
  
  r2 = r1 %>% mutate(SAMPLE_ID = row.names(r1)) %>%
    left_join(d %>% select(SAMPLE_ID, GENUS)) %>%
    mutate(genustype = ifelse(GENUS %in% c("Ctenotus", "Lerista"), GENUS, "other")) %>%
    group_by(genustype) %>%
    summarize(DR = mean(DR.t1.))
  
  rates = r2$DR
  names(rates) = r2$genustype
  
  ctdiff = rates["Ctenotus"] / rates["other"]
  lediff = rates["Lerista"] / rates["other"]
  
  res2[[i]] = gps %>% dplyr::select(-sp) %>%
    mutate(threshold = thresh[i])
  res[[i]] = c(Ntip(t1), thresh[i], gamma, ctdiff, lediff)
}

dd = data.frame(do.call("rbind", res))
names(dd) = c("Ntip", "threshold", "gamma", "ctdiff", "lediff")
dd1 = dd %>% select(threshold, ctdiff, lediff) %>%
  gather(type, value, -threshold)
x = data.frame(do.call("rbind", res2)) %>%
  spread(threshold, rates) %>%
  rename(SAMPLE_ID = ind) %>%
  left_join(d %>% select(SAMPLE_ID, OPERATIONAL_TAXON)) %>%
  group_by(OPERATIONAL_TAXON) %>% 
  sample_n(1) %>% ungroup()

c = ggplot(dd, aes(threshold, Ntip)) +
  geom_point() +
  geom_vline(xintercept = 2.5, col = "red", linetype = "dotted") +
  ylab("# of threshold species") +
  xlab(expression(tau ~ "(myr)")) +
  theme_classic()
a = ggplot(dd, aes(threshold, gamma)) +
  geom_point() +
  geom_vline(xintercept = 2.5, col = "red", linetype = "dotted") +
  ylab(expression(gamma)) +
  xlab(expression(tau ~ "(myr)")) +
  theme_classic()
b = ggplot(dd1, aes(threshold, value, color = type)) +
  geom_point(shape = 16) +
  geom_vline(xintercept = 2.5, col = "red", linetype = "dotted") +
  xlab(expression(tau ~ "(myr)")) +
  ylab(expression("relative difference in" ~ lambda[dr])) +
  scale_color_manual(values = c("black", "gray80"),
                     labels = c("Ctenotus",
                                "Lerista")) +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position = c(.8,.8))

c = ggplot(x, aes(`0.5`, `2.6`)) + 
  geom_point() + 
  xlab(expression(lambda[DR] ~ ", threshold 0.5 myr")) +
  ylab(expression(lambda[DR] ~ ", threshold 2.6 myr")) +
  theme_classic()
d = ggplot(x, aes(`2.6`, `5`)) + 
  geom_point() +
  xlab(expression(lambda[DR] ~ ", threshold 2.6 myr")) +
  ylab(expression(lambda[DR] ~ ", threshold 5 myr")) +
  theme_classic()

ab = plot_grid(plotlist = list(c, a, b), labels = c("A", "B", "C"), nrow = 1)
save_plot("figures/threshold_diversification.png", ab,
          base_height = 3, base_width = 12)
