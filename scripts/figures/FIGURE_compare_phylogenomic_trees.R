rm(list = ls())
library(ape)
library(phangorn)
library(phytools)

root_tree <- function(t, outs) {
  out2 = outs[outs %in% t$tip.label]
  t2 = root.phylo(t, outs[1])
  t2 = read.tree(text = write.tree(ladderize(t2)))
  tips = t2$tip.label
  tips = paste( d[match(tips, d$SAMPLE_ID), "SPECIES"], tips)
  tips[ match(out2, t2$tip.label) ] = 
    paste("outgroup", d[ match(out2, t2$tip.label), "SPECIES"])
  t2$tip.label = tips
  return(t2)
}

compare_tree <- function(tr1, tr2) {
 
  keep = intersect(tr1$tip.label, tr2$tip.label)
  tr1 = keep.tip(tr1, keep)
  tr2 = keep.tip(tr2, keep)
  
  nodes = Nnode(tr1)
  diffs = rep(NA, nodes)
  for (i in 1:nodes) {
    tips1 = tr1$tip.label[ Descendants(tr1, i + Ntip(tr1), type = "tips")[[1]] ]
    
    if (length(tips1) > 1) {
      focnode = findMRCA(tr2, tips1, type = "node")
      tips2 = tr2$tip.label[ Descendants(tr2, focnode, type = "tips")[[1]] ]
      diffs[i] = length(setdiff(tips2, tips1))
    } 
  }
  return(list(tr1, tr2, diffs))
}

process_tree <- function(tree) {
  aouts = outs[outs %in% tree$tip.label]
  tree1 = root(tree, aouts[1])
  tree2 = read.tree(text = write.tree(ladderize(drop.tip(tree1, aouts))))
  tree3 = drop.tip(tree2, drop_inds)
  tree3$tip.label[tree3$tip.label == "CCM5506"] = "NA_CCM5506_Ct_deca"
  tree3$tip.label[tree3$tip.label == "CCM6337"] = "NA_CCM6337_Ct_deca"
  
  keep = d %>% filter(d$SAMPLE_ID %in% tree3$tip.label) %>% group_by(OPERATIONAL_TAXON) %>%
    slice_sample() %>% dplyr::select(SAMPLE_ID) %>% pull()
  tree4 = keep.tip(tree3, keep)
  
  tree4$tip.label = d[match(tree4$tip.label, d$SAMPLE_ID), "OPERATIONAL_TAXON"]
  return(tree4)
}

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/phylogeny/")
d = read.csv("../metadata/sphenomorphine_all_v12.csv")
x = read.csv("../metadata/sphenomorphine_genomic_individuals_v1.csv")
outs1 = d[which(d$INGROUP == FALSE), "SAMPLE_ID"]
outs2 = x[which(x$INGROUP == FALSE), "SAMPLE_ID"]
outs = unique(c(outs1, outs2))

# two swapped inds
# why keep unknown samples?
drop_inds = c("CUMV_14452_Le_bipe", "UMMZ_244315_ct_quat", 
              "spheno_unk1", "spheno_unk3")

# load astral vs concat trees
at = read.tree("SqCL/astral/genetrees_miss0.3_tol1e-05_collapse80.tre")
ct = read.tree("SqCL/concat/concat_ind0.05_loci0.3_all_n5277.rooted.tre")

# root trees, ladderize, rename tips
t1 = process_tree(at)
t2 = process_tree(ct)

# compare astral & concat
res = compare_tree(t2, t1)
inv = compare_tree(t1, t2)

sum(res[[3]] > 0)

png("../../Sphenomorphine_Phylogeny2/figures/compare_phylogenomic.png",
    width = 15, height = 20, units = "in", res = 100)
par(xpd = T, mfrow = c(1, 2))
par(mar = c(0, 0, 0, 5))
# concat first

plot.phylo(res[[1]], show.tip.label = F)
title("concatenated", line = -2)
tiplabels(res[[1]]$tip.label, adj = c(0, 0.5),
          cex = 0.7, frame = "none")
for (i in 1:length(res[[3]])) {
  if (res[[3]][i] > 0) {
    nodelabels("", Ntip(res[[1]]) + i, pch = 16,
               col = "red", cex = 1, frame = "none")
  }
}
par(mar = c(0, 0, 0, 10))
plot.phylo(res[[2]], show.tip.label = F)
title("coalescent-based", line = -2)
tiplabels(res[[2]]$tip.label, adj = c(0, 0.5),
          cex = 0.7, frame = "none")
for (i in 1:length(inv[[3]])) {
  if (inv[[3]][i] > 0) {
    nodelabels("", Ntip(inv[[1]]) + i, pch = 16,
               col = "red", size = 1, frame = "none")
  }
}
dev.off()
