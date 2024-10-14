rm(list = ls())
library(ape)

# https://github.com/macroevolution/squamata/blob/main/scripts/6.supplemental_figures/SupplInfo_FigS6.R
# function that returns coordinates that outline a clade in a phylogeny
getCladePoly <- function(tree, node) {
  
  spanningInd <- range(which(tree$tip.label %in% geiger::tips(tree, node)))
  # spanningTips <- getSpanningTips(tree, node)
  
  # spanningInd <- which(tree$tip.label %in% spanningTips)
  
  nn <- min(spanningInd)
  
  nodePath <- nn
  while (nn != node) {
    nn <- tree$edge[tree$edge[,2] == nn, 1]
    nodePath <- c(nodePath, nn)
  }
  
  xx1 <- c()
  yy1 <- c()
  for (i in 1:length(nodePath)) {
    
    if (i > 1) {
      xx1 <- c(xx1, get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[nodePath[i]])
      desc <- phangorn::Descendants(tree, nodePath[i], type = 'children')
      yy1 <- c(yy1, min(get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[desc]))
    }
    
    xx1 <- c(xx1, get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[nodePath[i]])
    yy1 <- c(yy1, get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[nodePath[i]])
    
  }
  
  nn <- max(spanningInd)
  
  nodePath <- nn
  while (nn != node) {
    nn <- tree$edge[tree$edge[,2] == nn, 1]
    nodePath <- c(nodePath, nn)
  }
  
  xx2 <- c()
  yy2 <- c()
  for (i in 1:length(nodePath)) {
    
    if (i > 1) {
      xx2 <- c(xx2, get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[nodePath[i]])
      desc <- phangorn::Descendants(tree, nodePath[i], type = 'children')
      yy2 <- c(yy2, max(get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[desc]))
    }
    
    xx2 <- c(xx2, get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[nodePath[i]])
    yy2 <- c(yy2, get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[nodePath[i]])
    
  }
  
  xx <- c(xx1, rev(xx2), xx1[1])
  yy <- c(yy1, rev(yy2), yy1[1])
  
  return(cbind(xx, yy))
  
}

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation//")

t = read.tree("phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")
d = read.csv("metadata/sphenomorphine_all_v11.csv")
outs = d[which(d$INGROUP == FALSE), "SAMPLE_ID"]
t1 = drop.tip(t, outs)

# color constraints
edgecols = rep(alpha("gray50", 0.5), nrow(t1$edge))
x = read.tree("phylogeny/SqCL/concat/concat_ind0.05_loci0.3_all_n5277.phy.treefile.tre")
gen = x$tip.label[x$tip.label %in% t1$tip.label]
for (i in 1:length(gen)) {
  tip = which(t1$tip.label == gen[i])
  anc = c(tip, phangorn::Ancestors(t1, tip, type = "all"))
  edgecols[ which(t1$edge[, 2] %in% anc) ] = "dodgerblue"
}

pdf('../Sphenomorphine_Phylogeny2/figures/constraint.pdf', width = 13, height = 20)

par(mar = c(0,0,0,0), mfrow = c(1, 2))
plot.phylo(t1, show.tip.label = FALSE,
           edge.col = edgecols, cex = 0.5)

par(mar = c(0,0,0,5))
plot.phylo(t1, show.tip.label = FALSE, cex = 0.5, edge.color = alpha("gray50", 0.5))
pp = list.files('phylogeny/synthetic_trees//', 
                pattern = "all.constraint.treefile",
                full.names = T)
for (i in 1:length(pp)) {
  tips = read.tree(pp[i])$tip.label 
  tips = tips[tips %in% t1$tip.label]
  if (length(tips) > 1) {
  cladeMRCA <- getMRCA(t1, tips)
  coords <- getCladePoly(t1, cladeMRCA)
  
  polygon(coords, col = adjustcolor('red', alpha.f = 0.05), border = "red", lwd = 0.8)	
  
  spanningInd <- which(t1$tip.label %in% geiger::tips(t1, cladeMRCA))
  
  lab <- gsub(".*//", "", gsub(".all.constraint.treefile", "", pp[i]))
  text(x = max(branching.times(t1)), y = mean(spanningInd), lab, pos = 4, xpd = NA, cex = 0.75)
  }
}
dev.off()