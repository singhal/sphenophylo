rm(list = ls())
library(ape)

###########
# this creates a cluster level tree
###########

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")

x = read.csv("../Sphenomorphine_speciation/metadata/sphenomorphine_all_v12.csv")
t = read.tree("../Sphenomorphine_speciation/phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")
t = drop.tip(t, x[x$INGROUP == FALSE, "SAMPLE_ID"])
d = read.csv("data/structure/inds_to_clusters.29June24.csv")

## deal with taxonomic changes first
d$OPERATIONAL_TAXON = d$CURRENT_TAXON
d[d$OPERATIONAL_TAXON == "Concinnia_tenuis_1", "OPERATIONAL_TAXON"] = "Concinnia_tenuis"
d[d$OPERATIONAL_TAXON == "Lerista_microtis", "OPERATIONAL_TAXON"] = "Lerista_microtis_1"
d[d$OPERATIONAL_TAXON == "Lerista_zietzi", "OPERATIONAL_TAXON"] = "Lerista_chalybura"
d[d$OPERATIONAL_TAXON == "Notoscincus_wotjulum", "OPERATIONAL_TAXON"] = "Notoscincus_ornatus_4"
d$pop = gsub("Concinnia_tenuis_1", "Concinnia_tenuis", d$pop)
d$pop = gsub("Lerista_microtis", "Lerista_microtis_1", d$pop)
d$pop = gsub("Lerista_zietzi", "Lerista_chalybura", d$pop)
d$pop = gsub("Notoscincus_wotjulum", "Notoscincus_ornatus_4", d$pop)

# first make subsampled tree
inds1 = d %>% group_by(pop) %>%
  summarize(ind = sample(ind, 1)) %>% 
  pull(ind)

# add in individuals not in clusters
sps1 = unique(d$OPERATIONAL_TAXON)
sps1[!(sps1 %in% x$OPERATIONAL_TAXON)]
inds2 = x[!x$OPERATIONAL_TAXON %in% d$OPERATIONAL_TAXON, ] %>% 
  filter(complete.cases(OPERATIONAL_TAXON)) %>% 
  group_by(OPERATIONAL_TAXON) %>%
  summarize(SAMPLE_ID = sample(SAMPLE_ID, 1)) %>% 
  pull(SAMPLE_ID)
td = keep.tip(t, c(inds1, inds2))
td$tip.label = coalesce(d[match(td$tip.label, d$ind), "pop"], x[match(td$tip.label, x$SAMPLE_ID), "OPERATIONAL_TAXON"])
write.tree(td, "data/diversification/species_level_phylogeny.CLUSTERS.tre")

d2 = read.csv("data/operational_taxonomy/operational_taxonomy.csv")
d2$species = d2$OTU
d2$genus = gsub("_\\S+", "", d2$species)
tt = read.tree("data/diversification/species_level_phylogeny.OPERATIONAL_TAXON.tre")
d2$sampled = ifelse(d2$species %in% tt$tip.label, TRUE, FALSE)

cts = d2 %>% group_by(genus) %>%
  summarize(n = n()) %>% ungroup()
cts = d2 %>% filter(sampled == TRUE) %>% group_by(genus) %>%
  summarize(sampled = n()) %>% ungroup() %>% right_join(cts) %>%
  mutate(sampling_frac = round(sampled / n, 3))
dt = d2 %>% select(-sampled) %>% left_join(cts)

names = gsub("_V\\d+", "", td$tip.label)
names[!names %in% dt$species]

svec = dt[match(names, dt$species), "sampling_frac"]
table(complete.cases(svec))
out = paste0("data/diversification/", "sampling_fraction.CLUSTERS.csv")
write(svec, sep = "\n", file = out)
out2 = paste0("data/diversification/", "sampling_fraction.CLUSTERS.bamm.csv")
cat("1.0\n", file = out2)

x2 = data.frame(TIP = td$tip.label, TREE = names) %>%
  left_join(dt, by = c("TREE" = "species"))

write.table(x2 %>% select(TIP, genus, sampling_frac), out2,
              col.names = F, sep = "\t", row.names = F, quote = F,
              append = T)
