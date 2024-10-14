rm(list = ls())
library(VennDiagram)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")

t = read.tree("../Sphenomorphine_speciation/phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")
#######

d = read.csv("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/metadata/sphenomorphine_all_v10.csv")
x = list()
x[["mtDNA"]] = d[which(d$CYTB == TRUE), "SAMPLE_ID"]
x[["target-capture"]] = d[which(d$SQCL == TRUE), "SAMPLE_ID"]
x[["ddRAD"]] = d[which(d$ddRAD == TRUE & d$RECENT_RUN == FALSE), "SAMPLE_ID"]

for (type in names(x)) {
  tips = x[[type]]
  x[[type]] = tips[tips %in% t$tip.label]
}

length(unique(unlist(x)))

venn.diagram(x, filename = "figures/marker_types.png")
