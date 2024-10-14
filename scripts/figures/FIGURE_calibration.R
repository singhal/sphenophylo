rm(list = ls())

library(ape)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")

d = read.csv("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/metadata/sphenomorphine_all_v13.csv")
t = read.tree("../Sphenomorphine_speciation/phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")

d$otu_genus = gsub("_\\S+", "", d$OPERATIONAL_TAXON)
d1 = d %>% filter(complete.cases(otu_genus)) %>% group_by(otu_genus) %>% 
  slice_sample(n = 1) %>% select(SAMPLE_ID, otu_genus)
t1 = keep.tip(t, c(d1$SAMPLE_ID, d[which(d$INGROUP == FALSE), "SAMPLE_ID"]))
t1 = read.tree(text = write.tree(ladderize(t1)))

tipnames = coalesce(pull(d1[match(t1$tip.label, d1$SAMPLE_ID), "otu_genus"]), 
         d[match(t1$tip.label, d$SAMPLE_ID), "MORPHOLOGICAL_TAXON"])

par(mar = c(0, 0, 0, 12), xpd = T)
plot(t1, show.tip.label = F)
tiplabels(gsub("_", " ", tipnames), adj = c(0, 0.5), frame = "none", font = 3)

# OZ_SPHENO_L SAMR_53973_Ct_saxa AMSR_121997_An_swan
n1 = phytools::findMRCA(t1, tips = c("AMSR_121997_An_swan", "WAMR_141391_Ct_pant"))
nodelabels("", n1, frame = "none", pch = 21, cex = 2)
nodelabels("1", n1, frame = "none", cex = 1)

# SAHUL_SPHENO_J CAS_236743 AMSR_121997_An_swan
n2 = phytools::findMRCA(t1, tips = c("AMSR_121997_An_swan", "CAS_236743"))
nodelabels("", n2, frame = "none", pch = 21, cex = 2)
nodelabels("2", n2, frame = "none", cex = 1)

# ALLSPHENO_H CAS_213330 AMSR_121997_An_swan
n3 = phytools::findMRCA(t1, tips = c("AMSR_121997_An_swan", "CAS_213330"))
nodelabels("", n3, frame = "none", pch = 21, cex = 2)
nodelabels("3", n3, frame = "none", cex = 1)
