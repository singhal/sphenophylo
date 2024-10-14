rm(list = ls())

# script borrows (i.e., copies) heavily
# from Dan Rabosky's scripts
# https://github.com/macroevolution/squamata/blob/main/scripts/2.time_calibration_imputation/missing-tips-fxns.R

library(ape)
library(TreeSim)
library(phytools)

# Get elements of x that are not in y,
#   when x and y are numeric vectors.
#   Things are "same" when abs(x - y) < tol

setdiff_numeric <- function(x, y, tol = 1e-07){
  
  noMatch <- logical(length(x))
  for (ii in 1:length(x)){
    mm <- min(abs(x[ii] - y))
    if (mm > tol){
      noMatch[ii] <- TRUE		 
    }
  }
  return(x[noMatch])
}

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")

# add in a certain number of tips
nmiss = seq(20, 500, by = 20)
# do each simulation 5 times
nrep = 5

res = vector("list", length(nmiss) * nrep)

ct = 1
for (j in 1:length(nmiss)) {
  message("doing ", nmiss[j])
  for (x in 1:nrep) {
    # take tree
    phy = read.tree("data/diversification/species_level_phylogeny.OPERATIONAL_TAXON.tre")
    
    # max age of a new tip
    hts = data.frame(nodeHeights(phy)[which(phy$edge[,2] %in% 1:Ntip(phy)), ])
    hts$age = hts$X2 - hts$X1
    maxage = median(hts$age)
    
    # will use the global lambda to place tips
    lambda <- bd.ms(time = max(branching.times(phy)), n = Ntip(phy) + nmiss[j])
    
    # name the missing tips
    tipnames <- paste("baby", 1:nmiss[j], sep="")
    
    times = as.numeric(branching.times(phy))
    csbt  = corsim(times, lambda = lambda, mu=0, 
                   missing = nmiss[j], told=maxage)
    tree_age <- as.numeric(max(branching.times(phy)))
    stimes = tree_age - setdiff_numeric(csbt, times)
    
    for (ii in 1:length(tipnames)){
      phy <- BAMMtools:::getStartStopTimes(phy)
      goodset   <- phy$edge[,2][phy$begin <= stimes[ii] & phy$end > stimes[ii]]
      
      focal_node <- NA
      if (length(goodset) == 1){
        focal_node <- goodset
      }else{
        focal_node <- sample(goodset, 1)
      }
      
      edge_length = tree_age - stimes[ii]
      
      pos <- NA
      
      btimes = branching.times(phy)
      names(btimes) = seq(Ntip(phy) + 1, Ntip(phy) + Nnode(phy))
      
      if (focal_node <= length(phy$tip.label)){
        pos  <- edge_length
      }else{
        pos <- as.numeric(edge_length - btimes[as.character(focal_node)])
      }
      
      phy <- phytools:::bind.tip(phy, tip.label=tipnames[ii], where=focal_node, 
                                 edge.length=edge_length, position=pos)
    }
    res[[ct]] = c(nmiss[j], Ntip(phy), gammaStat(phy))
    ct = ct + 1
  }
}

res2 = data.frame(do.call("rbind", res))
write.csv(res2, "data/gamma/gamma_incipient_sim.csv")
