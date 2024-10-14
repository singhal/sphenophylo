rm(list = ls())
library(dplyr)
library(ape)
library(reshape2)

root_tree <- function(tt, sqcl) {
  outs = tt$tip.label[grep(".outgroup", tt$tip.label)]
  outs2 = gsub(".outgroup", "", outs)
  keep = sqcl$tip.label[ sqcl$tip.label %in% tt$tip.label ]
  outnodes = unlist(lapply(outs2, function(x) phytools::findMRCA(sqcl, c(x, keep))))
  outnodes2 = node.depth(sqcl)[outnodes]
  
  mainout = outs[ which(outnodes2 == max(outnodes2)) ]
  
  tt2 = root(tt, mainout[1], resolve.root = F)
  if (Ntip(tt2) > 1) {
    tt2 = read.tree(text = write.tree(ladderize(tt2)))
  }
  
  tt3 = drop.tip(tt2, outs)
  
  # get stem length
  if (Ntip(tt3) > 1) {
    stemnode = phytools::findMRCA(tt2, tt2$tip.label[grep("outgroup", tt2$tip.label, invert = T)])
    rootedge = tt2$edge.length[ which(tt2$edge[,2] == stemnode) ]
  } else {
    rootedge = tt2$edge.length[ which(tt2$edge[,2] == which(tt2$tip.label == tt3$tip.label)) ]
  }
  
  return(list(tt3, rootedge))
}

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

get_phy_dist <- function(tt) {
  m = cophenetic.phylo(tt[[1]])
  m = makeSymm(m)
  m1 = melt(as.matrix(m), varnames = c("ind1", "ind2"))
  m1 = m1 %>% rename(ddrad = value)
  
  m = cophenetic.phylo(tt[[2]])
  m = makeSymm(m)
  m2 = melt(as.matrix(m), varnames = c("ind1", "ind2"))
  m2 = m2 %>% rename(sqcl = value)
  
  mm = inner_join(m1, m2)
  # mm = mm[mm$ind1 != mm$ind2, ]
  return(mm)
}

rescale_trees <- function(t1, t2) {
  t1y = t1
  t1y$tip.label = gsub(".outgroup", "", t1y$tip.label)
  dd = get_phy_dist(list(t1y, t2))
  lm1 = lm(dd$sqcl ~ dd$ddrad)
  r2 = summary(lm1)$r.squared 
  sc = coef(lm1)[2]
  intercept = coef(lm1)[1]
  t1x = t1
  t1x$edge.length = t1$edge.length * sc
  # dd2 = get_phy_dist(list(t1x, t2))
  return(list(t1x, sc, r2, intercept))
}

process_sqcl <- function(sqcl) {
  sqcl$tip.label[which("CCM5506" == sqcl$tip.label)] = "NA_CCM5506_Ct_deca"
  sqcl$tip.label[which("CCM6337" == sqcl$tip.label)] = "NA_CCM6337_Ct_deca"
  drop = c("UMMZ_244315_ct_quat", "CUMV_14452_Le_bipe")
  sqcl = drop.tip(sqcl, drop)
  sqcl = drop.tip(sqcl, c("spheno_unk1", "spheno_unk3"))
  date = c("UMMZ_190536", "MVZ_Herp_256140",
           "UMMZ_242506", "CAS_236743",
           "AMNH_R153698", "UMFS_11369",
           "UMMZ_202015", "CAS_213330")
  outs = d[which(d$INGROUP == FALSE), "SAMPLE_ID"]
  sqcl = drop.tip(sqcl, outs[!outs %in% date])
  sqcl = root(sqcl, "CAS_213330", resolve.root = T)
  return(sqcl)
}

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation//")
d = read.csv("metadata/sphenomorphine_all_v11.csv")

# these individuals have ddRAD data but very low quality
dropinds = c("SAMR_56616_Ct_atla", "WAMR_135958_Ct_pant",
             "WAMR_98141_Le_pete", "SAMR_57650_Ct_atla")
# all outgroups are monophyletic
# except for C. burb -- what to do?
# tried to move it around, but not clear 
# where it belongs - just drop
burb = d[which(d$CURRENT_TAXON == "Ctenotus_burbidgei"), "SAMPLE_ID"]
dropinds = c(dropinds, burb)

type = "all"
sqcl = process_sqcl(read.tree("phylogeny/SqCL/concat/concat_ind0.05_loci0.3_all_n5277.phy.treefile.tre"))
t = sqcl 

pp = list.files('../Sphenomorphine_Phylogeny2/data/phylogenetic_support/', 
    pattern = "all.constraint.treefile",
    full.names = T)
    
for (i in 1:length(pp)) {
    x = read.tree(pp[i])
      
    # get rid of bad quality individuals
    x = drop.tip(x, dropinds)
    # rescale tree to match sqcl tree
    xx = rescale_trees(x, sqcl)
    x1 = xx[[1]]
    # root tree
    x2 = root_tree(x1, sqcl)
    if (length(x2[[2]]) == 0) {
        x2[[2]] = 0.001
    }
    root_edge = x2[[2]]
    x1 = x2[[1]]
        
    shared = x1$tip.label[ x1$tip.label %in% t$tip.label ]
    if (length(shared) > 1) {
      keep = sample(shared, 1)
      drop = shared[! shared %in% keep]
      t = drop.tip(t, drop)
    } else {
      keep = shared
    }
    node = phangorn::Ancestors(t, which(t$tip.label == keep))[1]
    t$tip.label[which(keep == t$tip.label)] = "new"
    x1$root.edge = root_edge
    t = bind.tree(t, x1, where = node)
    t = drop.tip(t, "new")
}

t1 = drop.tip(t, t$tip.label[!t$tip.label %in% d$SAMPLE_ID])
table(complete.cases(t1$node.label))
write.tree(t1, "~/Desktop/synthetic_tree.all.support.tree")