rm(list = ls())
library(tidyverse)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
d = read.csv("metadata/spheno_genera.csv")
x = readxl::read_xlsx("metadata/reptile_checklist_2024_03.xlsx")
x$genera = gsub(" \\S+", "", x$Species)
x1 = x[x$genera %in% d$genera, ]
x1$Species = gsub(" ", "_", x1$Species)

# now get in the changes that need to be recognized
d = readxl::read_xlsx("data/operational_taxonomy/sphenomorphine_genetic_papers.xlsx")
# confirm all my taxa are in there
pull(d[d$changes == "DROP", "taxon"]) %in% x1$Species
# then drop them
x2 = x1[! x1$Species %in% pull(d[d$changes == "DROP", "taxon"]), ]
# no one to add, so no worries there
d[d$changes == "ADD", ]

# now make sure all my operational units, as defined
# are in this operational taxonomy
x = read.csv("../Sphenomorphine_speciation/metadata/sphenomorphine_all_v12.csv")
# scrub off endings
otus = unlist(lapply(strsplit(x$OPERATIONAL_TAXON, "_"), function(x) paste0(x[1], "_", x[2])))
unique(otus[!(otus %in% x2$Species)])

# now print out the taxonomy
otu = data.frame(OTU = x$OPERATIONAL_TAXON, scrubbed = otus) %>% 
  distinct() %>% filter(complete.cases(OTU))

final = data.frame(OTU = c(otu$OTU, pull(x2[!(x2$Species %in% otu$scrubbed), "Species"])))
write.csv(final %>% arrange(OTU), "data/operational_taxonomy/operational_taxonomy.csv",
          row.names = F)