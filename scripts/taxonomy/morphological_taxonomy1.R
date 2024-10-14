rm(list = ls())
library(squamata)
library(tidyverse)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
d = read.csv("metadata/spheno_genera.csv")
x = readxl::read_xlsx("metadata/reptile_checklist_2024_03.xlsx")
x$genera = gsub(" \\S+", "", x$Species)
x1 = x[x$genera %in% d$genera, ]

d = readxl::read_xlsx("~/Desktop/ASH+Official+Species+List+30+November+2023+Simplified+Version+Without+Comments.xlsx")
x1[!x1$Species %in% d$species, ]
d[!d$species %in% x1$Species, "species"]


# start with current species, drop those that are "new"
# double check none were elevated from subspecies
###  to make a list of species that are only morphological

###  confirm no cases of synonymization
syn = vector("list", length = nrow(x1))
for (i in 1:nrow(x1)) {
  sp = pull(x1[i, "Species"])
  dd = squamata::reptileDB_synonymsFromAccepted(sp)[[1]]
  dd$species = sp
  
  num_recent = dd %>% filter(year > 1960) %>% nrow()
  if (num_recent > 0) {
    dd = dd %>% filter(year > 1960)
  }
  
  dd$nonsyn = ''
  num_syn = length(unique(dd$synonym))
  if (num_syn == 1) {
    dd$nonsyn = TRUE
  } 
  
  syn[[i]] = dd
}
syn2 = do.call("rbind", syn)
write.csv(syn2, "~/Desktop/all_synonym_sphenomorphines.csv", row.names = F)

# to possibly identify synonymized species
# get all synonyms
syns = gsub("_", " ", unique(syn2$synonym))
syn = data.frame(genus = unlist(lapply(strsplit(syns, " "), function(x) x[[1]])),
                 sp = unlist(lapply(strsplit(syns, " "), function(x) x[[2]])))
syn = syn %>% filter(genus %in% d$genera) %>% distinct() %>%
  mutate(fullname = paste(genus, sp)) %>% 
  filter(!fullname %in% x1$Species)
write.csv(syn, "~/Desktop/synonym_sphenomorphines.csv", row.names = F)

x1$year = parse_number(x1$Author)
write.csv(x1, "~/Desktop/starting_morphological_taxonomy.csv", row.names = F)
