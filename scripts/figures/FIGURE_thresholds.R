rm(list = ls())
library(ape)
library(phytools)
library(PCPS)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation//")
t = read.tree("phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")
d = read.csv("metadata/sphenomorphine_all_v13.csv")

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

tt = unique(d[complete.cases(d$CURRENT_TAXON), "CURRENT_TAXON"])
thresh = seq(0.1, 10, 0.1)
res = vector("list", length(thresh))
for (i in 1:length(thresh)) {
  sp = get_groups(t, thresh[i])
  
  match = data.frame(lineage = tt,
                     num_sps1 = NA,
                     num_inds1 = NA,
                     num_sps2 = NA,
                     num_inds2 = NA,
                     threshold = NA)
  for (j in 1:length(tt)) {
    sp1 = tt[j]
    inds1 = d[which(d$CURRENT_TAXON == sp1), "SAMPLE_ID"]
    sp2 = unique(sp[sp$ind %in% inds1, "sp"])
    inds2 = sp[sp$sp %in% sp2, "ind"]
    sp1 = unique(d[which(d$SAMPLE_ID %in% inds2), "CURRENT_TAXON"])
    
    match[j, ] = c(tt[j], length(sp1), length(inds1), length(sp2), length(inds2), thresh[i])
  }
  res[[i]] = match
  cat(thresh[i], "\n")
}

res2 = do.call("rbind", res)
res2$same = TRUE
res2[res2$num_sps1 != res2$num_sps2, "same"] = FALSE
res2[res2$num_inds1 != res2$num_inds2, "same"] = FALSE
res3 = res2 %>% group_by(threshold) %>% 
  summarize(percent_same = sum(same) / length(same)) %>% 
  ungroup() %>% mutate(threshold = as.numeric(threshold))

a = ggplot(res3, aes(threshold, percent_same)) + geom_point() +
  ylab("% species identical") +
  xlab("phylogenetic threshold (Myr)") +
  theme_classic()
cowplot::save_plot("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/figures/threshold.png",
                   a, base_height = 4, base_width = 6)
