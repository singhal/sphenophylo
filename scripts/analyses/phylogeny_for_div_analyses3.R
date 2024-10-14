rm(list = ls())

library(ape)
library(phytools)
library(tidyverse)

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

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")

d = read.csv("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/metadata/sphenomorphine_all_v13.csv")

t = read.tree("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")
t = drop.tip(t, d[d$INGROUP == FALSE, "SAMPLE_ID"])

d = d[d$SAMPLE_ID %in% t$tip.label, ]

xx = get_groups(t, 2.5) 
xx$GENUS = gsub("_\\S+", "", d[match(xx$ind, d$SAMPLE_ID), "OPERATIONAL_TAXON"])
xx2 = xx %>% group_by(sp, GENUS) %>% summarize(n = n()) %>% 
  slice_max(n) %>% ungroup() %>% mutate(SPECIES = paste0(GENUS, "_", sp))
xx = left_join(xx, xx2)

d$THRESHOLD_TAXON = xx[match(d$SAMPLE_ID, xx$ind), "SPECIES"]
d1 = d[complete.cases(d[, "THRESHOLD_TAXON"]), ]
t1 = keep.tip(t, d1$SAMPLE_ID)
tips = t1$tip.label
  
sps = d1[match(tips, d1$SAMPLE_ID), "THRESHOLD_TAXON"]
names(sps) = tips
usps = unique(sps)
for (i in 1:length(usps)) {
    sp = usps[i]
    inds = names(sps[ which(sps == sp) ])
    if (length(inds) > 1) {
      desc = t1$tip.label[ phangorn::Descendants(t1, phytools::findMRCA(t1, inds), type = "tips")[[1]] ]
      ndiff = setdiff(desc, inds)
      if (length(ndiff) > 0) {
        cat(type, sp, "\n")
        alldesc = phytools::getDescendants(t1, phytools::findMRCA(t1, inds))
        alldesc2 = lapply(alldesc, function(x) t1$tip.label[ phangorn::Descendants(t1, x, type = "tips")[[1]] ])
        alldesc3 = alldesc2[ which(unlist(lapply(alldesc2, function(x) length(setdiff(x, inds)))) == 0) ]
        cts = unlist(lapply(alldesc3, function(x) length(x)))
        ind = sample(alldesc3[[ which(cts == max(cts))[1] ]], 1)
      } else {
        ind = sample(inds, 1)
      }
      drop = inds[!inds %in% ind] 
      t1 = drop.tip(t1, drop)
    } 
  }
t1$tip.label = sps[ t1$tip.label ]
out = paste0("data/diversification/", "species_level_phylogeny.THRESHOLD_TAXON.tre")
write.tree(t1, out)

# sampling frac
# use operational taxonomy because best we can do
d2 = read.csv("data/operational_taxonomy/operational_taxonomy.csv")
d3 = d2 %>% group_by(genus) %>% summarize(n = length(species)) %>% ungroup()
d4 = xx2 %>% group_by(genus = GENUS) %>% summarize(sampled = length(sp)) %>% ungroup()
d5 = d3 %>% left_join(d4) %>% 
  mutate(sampled = ifelse(complete.cases(sampled), sampled, 0)) %>%
  mutate(sampling_frac = round(sampled / n, 3)) %>%
  mutate(sampling_frac = ifelse(sampling_frac > 1, 1, sampling_frac))
dt = xx2 %>% left_join(d5 %>% select(GENUS = genus, sampling_frac))


svec = pull(dt[match(t1$tip.label, dt$SPECIES), "sampling_frac"])
cat(table(complete.cases(svec)))
out = paste0("data/diversification/", "sampling_fraction.THRESHOLD_TAXON.csv")
write(svec, sep = "\n", file = out)
out2 = paste0("data/diversification/", "sampling_fraction.THRESHOLD_TAXON.bamm.csv")
cat("1.0\n", file = out2)
write.table(dt %>% filter(SPECIES %in% t1$tip.label) %>% select(SPECIES, GENUS, sampling_frac), out2,
              col.names = F, sep = "\t", row.names = F, quote = F,
              append = T)
