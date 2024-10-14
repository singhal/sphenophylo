rm(list = ls())
library(coda)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")

mcmcout <- read.csv("data/diversification/MORPHOLOGICAL_TAXON_mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
