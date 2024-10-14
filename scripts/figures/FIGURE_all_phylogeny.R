rm(list = ls())
library(sf)
library(ggplot2)
library(tidyverse)
library(ape)
library(phytools)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
source("scripts/colors.R")

#######
#	called by colorClade(...)
#
zip = function(x,y) {
  ret = numeric(2*length(x));
  i = 1; j =1;
  while(i <= length(ret)) {
    ret[i] = x[j];
    i = i+1;
    ret[i] = y[j];
    i = i+1;
    j = j+1;
  }
  return(ret);
}

#######
#	called by colorClade(...)
#
zip = function(x,y) {
  ret = numeric(2*length(x));
  i = 1; j =1;
  while(i <= length(ret)) {
    ret[i] = x[j];
    i = i+1;
    ret[i] = y[j];
    i = i+1;
    j = j+1;
  }
  return(ret);
}

########
# polar tree plotting function that returns coordinates of edges
#
#
polar.phylo = function(phy, tip_names, vtheta = 1, rbf = 0.001, labels = FALSE, lwd = 1, edge.color = 1, arc.color = 1, cex = 1)
{
  arc.color = edge.color
  phy = BAMMtools:::getStartStopTimes(phy)
  tH = max(phy$end);
  rb = tH * rbf;
  ret = BAMMtools:::setPolarTreeCoords(phy, vtheta, rbf);
  
  x0 = ret$segs[, 1];
  x1 = ret$segs[, 3];
  y0 = ret$segs[, 2];
  y1 = ret$segs[, 4];
  ofs = 0;
  if (labels) {
    ofs = max(nchar(tip_names) * 0.03 * cex);
  }
  
  plot.new();
  plot.window(xlim = c(-1, 1) + c(-rb, rb) + c(-ofs, ofs), ylim = c(-1, 1) + c(-rb, rb) + c(-ofs, ofs), asp = 1);
  segments(x0, y0, x1, y1, col = edge.color, lwd = lwd, lend = 2);
  BAMMtools:::arc(0, 0, ret$arcs[, 1], ret$arcs[, 2], c(rb, rb + phy$end/tH), border = rep(arc.color, nrow(ret$arcs)), lwd = lwd);
  
  invisible(ret);
}

#######
#	Arguments:
#		phy     -  phylo object 
#		p		-  plotting coordinates of phylogeny. returned invisibly from plot.bammdata
#		node1   -  left bound tip node (may be character string) 
#		node2   -  right bound tip node (may be character string)
#		col     -  fill and border color of polygon
#		alpha   -  transparency of color
colorClade = function(phy, p, nodes, col, border, alpha, name, cex) {
  n2 = max(nodes);
  n1 = min(nodes);
  seq.nod = .Call("seq_root2tip",phy$edge,phy$Nnode+1,phy$Nnode,PACKAGE="ape");
  nmrca = BAMMtools:::getmrca(phy, n1, n2);
  left = seq.nod[[n2]];
  left = left[which(left == nmrca):length(left)];
  right = seq.nod[[n1]];
  right = right[which(right == nmrca):length(right)];
  if (n1 == n2) {
    tips = seq.int(n1, n2-1);	
  }
  else if (n1 == n2+1) {
    tips = seq.int(n1, n2-1);
  }
  else {
    tips = seq.int(n1, n2-1);
  }
  
  lc = p[rownames(p) %in% left,];
  rc = p[rownames(p) %in% right,];
  tc = p[rownames(p) %in% tips,];
  
  if (is.null(dim(tc))) {
    tc = matrix(tc, nrow = 1);
  }
  
  xv = c(zip(lc[-1,1],lc[-1,3]), tc[nrow(tc):1,3], zip(rc[nrow(rc):2,3], rc[nrow(rc):2,1]));
  yv = c(zip(lc[-1,2],lc[-1,4]), tc[nrow(tc):1,4], zip(rc[nrow(rc):2,4], rc[nrow(rc):2,2]));
  polygon(xv, yv, col=BAMMtools:::transparentColor(col,alpha), border= border);
  
  angle = (mean(nodes) * 2 * 3.14159) / length(phy$tip.label)
  text(x=1.05 * cos(angle),
       y = 1.05 * sin(angle), name, srt=(angle * 360 / (2 * 3.14159)), 
       adj=c(0,0.5), cex=cex, font = 3)
}

d = read.csv("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/metadata/sphenomorphine_all_v10.csv")
t = read.tree("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")

outs = d[which(d$INGROUP == FALSE), "SAMPLE_ID"]
t1 = drop.tip(t, outs)
d1 = d[!d$SAMPLE_ID %in% outs, ]

colors2 = rep(allcols, 20)

pdf("figures/lineage_phylogeny.pdf", width=6, height=6)
par(mar=c(0, 0, 0, 0), xpd = TRUE)

cexval = 0.4
ret = polar.phylo(t1, lwd=0.5, vtheta=0, 
                  labels=FALSE, cex=cexval, edge.color = "gray60")

lins2 = unique(d[match(t1$tip.label, d$SAMPLE_ID), "CURRENT_TAXON"])
groups = split(d, d$CURRENT_TAXON)
for (i in 1:length(lins2)) {
  group = lins2[i]
  inds = groups[[group]]$SAMPLE_ID
  nodes = match(inds, t1$tip.label)
  group2 = gsub("_", " ", group)
  colorClade(t1, ret$seg, nodes, colors2[i], NA, 0.5, "", cexval)
}
dev.off()