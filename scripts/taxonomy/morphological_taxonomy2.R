rm(list = ls())

library(tidyverse)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")

d1 = read.csv("data/morphological_taxonomy/starting_morphological_taxonomy.csv")
d2 = read.csv("data/morphological_taxonomy/synonym_sphenomorphines.csv")

# check what I'm dropping
(drop = d1[d1$described_by_morphology == FALSE, ])
# drop them
dd = d1[d1$described_by_morphology == TRUE, ]

# add in any synonyms
d2[d2$synonym == TRUE, "fullname"]
# well, there aren't any

# save the new taxonomy
write.csv(dd %>% select(Species, Author, year),
          "data/morphological_taxonomy/morphological_taxonomy.retained.csv",
          row.names = FALSE)
# record what was dropped for sanity
write.csv(drop %>% select(Species, Author, year), 
          "data/morphological_taxonomy/morphological_taxonomy.dropped.csv",
          row.names = FALSE)
