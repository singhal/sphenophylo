rm(list = ls())
library(ape)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")

d = read.csv("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/metadata/sphenomorphine_all_v12.csv")

t = read.tree("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")
t = drop.tip(t, d[d$INGROUP == FALSE, "SAMPLE_ID"])

types = c("MORPHOLOGICAL_TAXON", "OPERATIONAL_TAXON")
res = vector("list", length(types))
names(res) = types

d = d[d$SAMPLE_ID %in% t$tip.label, ]

for (type in types) {

  d1 = d[complete.cases(d[, type]), ]
  t1 = keep.tip(t, d1$SAMPLE_ID)
  tips = t1$tip.label
  
  sps = d1[match(tips, d1$SAMPLE_ID), type]
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
  res[[type]] = t1
  out = paste0("data/diversification/", "species_level_phylogeny.", type, ".tre")
  write.tree(t1, out)
}
saveRDS(res, "data/diversification/lineage_level_phylogeny.Rds")

# sampling frac
genera = read.csv("metadata/spheno_genera.csv")
d1 = read.csv("data/morphological_taxonomy/morphological_taxonomy.retained.csv")
d1$species = gsub(" ", "_", d1$Species)
d1$genus = gsub("_\\S+", "", d1$species)
d1$type = "MORPHOLOGICAL_TAXON"
d1$sampled = ifelse(d1$species %in% res[[d1$type[1]]]$tip.label, TRUE, FALSE)

d2 = read.csv("data/operational_taxonomy/operational_taxonomy.csv")
d2$species = d2$OTU
d2$genus = gsub("_\\S+", "", d2$species)
d2$type = "OPERATIONAL_TAXON"
res[[d2$type[1]]]$tip.label %in% d2$species 
d2$sampled = ifelse(d2$species %in% res[[d2$type[1]]]$tip.label, TRUE, FALSE)

dd = rbind(d1 %>% select(species, genus, type, sampled),
      d2 %>% select(species, genus, type, sampled))

for (tt in types) {
  t = res[[tt]]
  cts = dd %>% filter(type == tt) %>% group_by(genus) %>%
    summarize(n = n()) %>% ungroup()
  cts = dd %>% filter(type == tt, sampled == TRUE) %>% group_by(genus) %>%
    summarize(sampled = n()) %>% ungroup() %>% right_join(cts) %>%
    mutate(sampled = ifelse(complete.cases(sampled), sampled, 0)) %>%
    mutate(sampling_frac = round(sampled / n, 3))
    
  dt = dd %>% filter(type == tt) %>% select(-sampled) %>% left_join(cts)

  svec = dt[match(t$tip.label, dt$species), "sampling_frac"]
  cat(table(complete.cases(svec)))
  out = paste0("data/diversification/", "sampling_fraction.", tt, ".csv")
  write(svec, sep = "\n", file = out)
  out2 = paste0("data/diversification/", "sampling_fraction.", tt, ".bamm.csv")
  cat("1.0\n", file = out2)
  write.table(dt %>% filter(species %in% t$tip.label) %>% select(species, genus, sampling_frac), out2,
              col.names = F, sep = "\t", row.names = F, quote = F,
              append = T)
}  
