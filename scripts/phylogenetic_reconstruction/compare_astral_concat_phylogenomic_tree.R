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

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/phylogeny/")
d = read.csv("../metadata/sphenomorphine_genomic_individuals_v1.csv")

outs = d[d$INGROUP == FALSE, "SAMPLE_ID"]

# load trees
t1 = list.files("SqCL/astral/", pattern = ".tre", full.names = T)
t2 = list.files("SqCL/concat/", pattern = ".treefile", full.names = T)
tt = c(t1, t2)
tt2 = lapply(tt, read.tree)

# root trees
tt3 = lapply(tt2, function(x) root_tree(x, outs))

# compare within astral - miss
res1 = compare_tree(tt3[[2]], tt3[[6]])
res2 = compare_tree(tt3[[4]], tt3[[6]])

# compare within astral - AHE
res3 = compare_tree(tt3[[5]], tt3[[6]])

# compare within concat - AHE
res4 = compare_tree(tt3[[7]], tt3[[8]])

# compare astral & concat
res5 = compare_tree(tt3[[6]], tt3[[8]])


pdf("~/Desktop/concatenated_all.pdf", width = 8, height = 20)
par(mar = c(0, 0, 0, 10), xpd = T)
t = tt3[[8]]
t2 = drop.tip(t, t$tip.label[ grep("outgroup", t$tip.label) ])
t2 = drop.tip(t2, "NA spheno_unk3")
plot.phylo(t2, show.tip.label = F)
tiplabels(t2$tip.label, frame = "none", adj = c(0, 0.5), cex = 0.5)
dev.off()

pdf("~/Desktop/ASTRAL_all.pdf", width = 8, height = 20)
par(mar = c(0, 0, 0, 10), xpd = T)
t = tt3[[6]]
t2 = drop.tip(t, t$tip.label[ grep("outgroup", t$tip.label) ])
t2 = drop.tip(t2, c("NA spheno_unk1", "NA spheno_unk3"))
plot.phylo(t2, show.tip.label = F)
tiplabels(t2$tip.label, frame = "none", adj = c(0, 0.5), cex = 0.5)
dev.off()

t1 = tt3[[8]]
t1 = drop.tip(t1, t1$tip.label[ grep("outgroup", t1$tip.label) ])
t2 = tt3[[6]]
t2 = drop.tip(t2, t2$tip.label[ grep("outgroup", t2$tip.label) ])

t1 = keep.tip(t1, intersect(t1$tip.label, t2$tip.label))
t2 = keep.tip(t2, intersect(t1$tip.label, t2$tip.label))

rel = cbind(t1$tip.label, t2$tip.label)
obj = cophylo(tr1=t1, tr2=t2)
pdf("~/Desktop/concatenated_astral_AHE_comparison.pdf", width=20, height=30)
plot(obj)
dev.off()
