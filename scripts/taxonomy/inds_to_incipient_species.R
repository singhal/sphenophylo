rm(list = ls()) 

library(sf)
library(tidyverse)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
d = list.files("data/structure/fastStructure/", pattern = "csv", full.names = T)
xx = read.csv("../Sphenomorphine_speciation/metadata/sphenomorphine_all_v11.csv")

res = vector("list", length(d))
for (i in 1:length(d)) {
  dd = read.csv(d[i])
  
  # use components to pick best K
  k = dd[2, 3]
  
  # get assignment values
  sp = gsub(".choose.*", "", gsub(".*//", "", d[i]))
  x = read.table(paste0("data/structure//fastStructure/", sp, ".", k, ".meanQ"))
  
  # ind names
  inds = read.table(paste0("data/structure/admixture/", sp, ".miss0.7.MAC2.fam"))
  x$ind = inds$V2
  
  res[[i]] = x %>% gather(key = pop, value = prob, -ind) %>%
    group_by(ind) %>% slice_max(prob) %>% 
    left_join(xx %>% select(ind = SAMPLE_ID, CURRENT_TAXON)) %>% 
    filter(complete.cases(CURRENT_TAXON)) %>% 
    mutate(pop = paste0(CURRENT_TAXON, "_", pop))
}

res2 = do.call("rbind", res)
write.csv(res2, "data/structure/inds_to_clusters.29June24.csv",
          quote = F, row.names = F)
