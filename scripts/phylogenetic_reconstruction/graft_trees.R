rm(list = ls())
library(dplyr)
library(ape)
library(RColorBrewer)
library(reshape2)

plot_tree <- function(t2, type, name) {
  sps = d[match(t2$tip.label, d$SAMPLE_ID), "CURRENT_TAXON"]
  tips = paste(sps,
               d[match(t2$tip.label, d$SAMPLE_ID), "CLUSTER_SUBGENUS"],
               t2$tip.label, sep = " / ")
  
  sps2 = unique(sps)
  colors = rep(NA, length(sps))
  for (j in 1:length(sps2)) {
    colors[ which(sps == sps2[j]) ] = allcols[j]
  }
  
  write.tree(t2, paste0("~/Desktop/synthetic_tree.", name, ".", type, ".rooted.18April2023.tre"))
  
  pdf(paste0("~/Desktop/synthetic_tree.", type, ".rooted.18April2023.pdf"), height = 360, width = 10)
  par(mar = c(0, 0, 0, 15), xpd = T)
  plot(t2, show.tip.label = F)
  tiplabels(tips, bg = colors, frame = "rect", adj = c(0, 0.5))
  dev.off()
}

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
d = read.csv("metadata/sphenomorphine_all_v8.csv")
allcols = rep(brewer.pal(12, "Set3"), 30)
types = c("all", "nuc")

# these individuals have ddRAD data but very low quality
dropinds = c("SAMR_56616_Ct_atla", "WAMR_135958_Ct_pant",
             "WAMR_98141_Le_pete", "SAMR_57650_Ct_atla")
# all outgroups are monophyletic
# except for C. burb -- what to do?
# tried to move it around, but not clear 
# where it belongs - just drop
burb = d[which(d$CURRENT_TAXON == "Ctenotus_burbidgei"), "SAMPLE_ID"]
dropinds = c(dropinds, burb)

sqcltrees = list(
  undated = process_sqcl(read.tree("phylogeny/SqCL/concat/concat_ind0.05_loci0.3_all_n5277.phy.treefile.tre")), 
  dated = process_sqcl(read.tree("phylogeny/SqCL/concat/concat_ind0.05_loci0.3_all_n5277.rooted.dated.tre")) 
)

# #######################
# # check that they match constraint tree
# #######################
# cc = data.frame(trees = pp, num_same = NA, nnode = NA)
# for (i in 1:length(pp)) {
#   x = read.tree(pp[i])
#   con = gsub("synthetic_trees", "constraints", pp[i])
#   con = gsub("\\.\\S+", ".tre", con)
#  
#   if (file.exists(con)) {
#     con = unroot(read.tree(con))
#     comp = unroot(keep.tip(x, con$tip.label))
#     
#     cout = ape::comparePhylo(con, comp)
#     cc[i, "num_same"] = as.numeric(gsub(" splits.*", "", cout$messages[7]))
#     cc[i, "nnode"] = Nnode(con)
#   }
# }
# cc$num_diff = cc$nnode - cc$num_same
# View(cc)
# #######################
# # great! by switching to iqtree v1
# # constraints were taken into account
# #######################
# 
# #######################
# # check that outgroups are mono
# #######################
# for (i in 1:length(pp)) {
#   x = read.tree(pp[i])
#   outs = x$tip.label[ grep("outgroup", x$tip.label)]
#   noouts = x$tip.label[ grep("outgroup", x$tip.label, invert = T)]
#   x1 = root(x, outs[1], resolve.root = TRUE)
#   if (length(noouts) > 1) {
#     mrca = phytools::findMRCA(x1, noouts, "node")
#     mrcaout = x1$tip.label[ phangorn::Descendants(x1, mrca, type = "tips")[[1]] ]
#     outs2 = outs[ outs %in% mrcaout ]
#     if (length(outs2) > 0) {
#       outx = gsub("phylogeny/synthetic_trees//", "~/Desktop/", pp[i])
#       outx = gsub("treefile", "tre", outx)
#       write.tree(x1, outx)
#     }
#   }
# }

res = vector("list", 4)
ix = 1
for (z in 1:length(sqcltrees)) {
  for (y in 1:length(types)) {
    type = types[y]
    sqcl = sqcltrees[[z]]
    t = sqcl
    name = names(sqcltrees)[z]
    
    pp = list.files('phylogeny/synthetic_trees/', 
                    pattern = paste0(type, ".constraint.treefile"),
                    full.names = T)
    
    rr = data.frame(tree = pp, rescale = NA, root_edge = NA, 
                    r2 = NA, type = type, name = name,
                    intercept = NA)
    
    for (i in 1:length(pp)) {
      x = read.tree(pp[i])
      
      # get rid of bad quality individuals
      x = drop.tip(x, dropinds)
      # rescale tree to match sqcl tree
      xx = rescale_trees(x, sqcl)
      x1 = xx[[1]]
      rr[i, "rescale"] = xx[[2]]
      rr[i, "r2"] = xx[[3]]
      rr[i, "intercept"] = xx[[4]]
      # root tree
      x2 = root_tree(x1, sqcl)
      if (length(x2[[2]]) == 0) {
        x2[[2]] = 0.001
      }
      root_edge = x2[[2]]
      rr[i, "root_edge"] = root_edge
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
    res[[ix]] = list(t, rr)
    ix = ix + 1
  }
}
# 
# # as a sanity check
# # look at how much root edge, r2 & rescale scalar
# # varies across trees
rr = do.call("rbind", lapply(res, function(x) x[[2]]))
rr$tree = gsub("\\..+", "", rr$tree)
rr$tree = gsub(".*//", "", rr$tree)
rr1 = rr %>% select(tree, type, name, r2) %>%
  pivot_wider(id_cols = tree,
              names_from = c("type", "name"),
              values_from = r2)
rr2 = rr %>% select(tree, type, name, root_edge) %>%
  pivot_wider(id_cols = tree,
              names_from = c("type", "name"),
              values_from = root_edge)
rr3 = rr %>% select(tree, type, name, rescale) %>%
  pivot_wider(id_cols = tree,
              names_from = c("type", "name"),
              values_from = rescale)
rr4 = rr %>% select(tree, type, name, intercept) %>%
  mutate(intercept = abs(intercept)) %>%
  pivot_wider(id_cols = tree,
              names_from = c("type", "name"),
              values_from = intercept)
View(rr4)

# # now for each tree
# # look at how much distances vary relative to 
# # original tree
# # and sqcl tree
# # measure as correlations 
# #   and residuals to 1 : 1 line
ix = 1
res2 = vector("list", 4)
for (z in 1:length(sqcltrees)) {
  for (y in 1:length(types)) {
    type = types[y]
    sqcl = sqcltrees[[z]]
    name = names(sqcltrees)[z]
    restree = res[[ix]][[1]]

    pp = list.files('phylogeny/synthetic_trees/',
                    pattern = paste0(type, ".constraint.treefile"),
                    full.names = T)

    rr = data.frame(tree = pp,
                    r2 = NA,
                    residual = NA,
                    intercept = NA,
                    slope = NA,
                    num = NA,
                    type = type,
                    name = name)
    pdf(paste0("~/Desktop/pairwise.", type, ".", name, ".pdf"),
        height = 24, width = 12)
    par(mfrow = c(8, 4), mar = c(3, 3, 3, 1))
    for (i in 1:length(pp)) {
      x = read.tree(pp[i])

      # get rid of bad quality individuals
      x = drop.tip(x, dropinds)

      x$tip.label = gsub(".outgroup", "", x$tip.label)
      shared = intersect(x$tip.label, sqcl$tip.label)

      tt1 = keep.tip(restree, shared)
      tt2 = keep.tip(sqcl, shared)

      pd = get_phy_dist(list(tt1, tt2))
      lm1 = lm(pd$ddrad ~ pd$sqcl)
      plot(pd$sqcl, pd$ddrad)
      abline(a = 0, b = 1, col = "red")

      pltix = gsub("\\..+", "", pp[i])
      pltix = gsub(".*//", "", pltix)
      mtext(pltix, side = 3, line = 1)

      rr[i, "residual"] = abs(sum(pd$ddrad - pd$sqcl) / length(shared))
      rr[i, "num"] = length(shared)
      rr[i, "r2"] = summary(lm1)$r.squared
      rr[i, "intercept"] = coef(lm1)[1]
      rr[i, "slope"] = coef(lm1)[2]
    }

    rr2 = get_phy_dist(list(sqcl, res[[ix]][[1]]))
    rr2$type = type
    rr2$name = name

    res2[[ix]] = list(rr, rr2)
    ix = ix + 1
    dev.off()
  }
}

rr2 = do.call("rbind", lapply(res2, function(x) x[[1]]))
rr2$tree = gsub("\\..+", "", rr2$tree)
rr2$tree = gsub(".*//", "", rr2$tree)
rr2a = rr2 %>% select(tree, type, name, r2) %>%
  pivot_wider(id_cols = tree,
              names_from = c("type", "name"),
              values_from = r2)
rr2b = rr2 %>% select(tree, type, name, slope) %>%
  pivot_wider(id_cols = tree,
              names_from = c("type", "name"),
              values_from = slope)
rr2c = rr2 %>% mutate(intercept = abs(intercept)) %>%
  select(tree, type, name, intercept) %>%
  pivot_wider(id_cols = tree,
              names_from = c("type", "name"),
              values_from = intercept)
View(rr2c)

rr3 = res2[[1]][[2]]
rr3 = rr3[rr3$ddrad > 0, ]
lm1 = lm(rr3$ddrad ~ rr3$sqcl)
rr3$predict_ddrad = predict(lm1, data = rr3$sqcl)
rr3$error = abs(rr3$predict_ddrad - rr3$ddrad) / rr3$predict_ddrad
hist(rr3$error, breaks = 1000)

################################## 
# root edge missing for inornatus
#   everyone else is fine
# big variance in rescale - 
##################################
# 
# ###############################
# # check monophyly of taxon
# ###############################
# sps = d[match(t$tip.label, d$SAMPLE_ID), "CURRENT_TAXON"]
# names(sps) = t$tip.label
# sps2 = unique(sps)
# for (i in 1:length(sps2)) {
#   inds = names(sps[which(sps == sps2[i])])
#   if (length(inds) > 1) {
#     node = phytools::findMRCA(t, inds)
#     inds2 = t$tip.label[ phangorn::Descendants(t, node, type = "tips")[[1]] ]
#     allsp = unique(sps[ inds2 ])
#     if (length(allsp) > 1) {
#       cat(sps2[i], paste(allsp, collapse = ": "), "\n")
#     }
#   }
# }
# ###############################
# # all taxon are monophyletic!
# ###############################
# 
# 
# ###############################
# # check long branch lengths
# ###############################
# t2 = read.tree(text = write.tree(ladderize(t)))
# edgewid = rep(0.1, length(t2$edge.length))
# edgecols = rep("black", length(t2$edge.length))
# edgecols[ which(t2$edge.length > 0.01) ] = "red"
# edgewid[ which(t2$edge.length > 0.01) ] = 2
# pdf(paste0("~/Desktop/branch_lengths.", type, ".pdf"), height = 360, width = 10)
# par(mar = c(0, 0, 0, 0), xpd = T)
# plot(t2, show.tip.label = T, edge.color = edgecols, edge.width = edgewid)
# dev.off()
# ###############################
# # most of the deep nodes are where you would expect
# # the two that are concerning - regius & notoscincus
# ###############################

# #######################
# # missing inds
# #######################
# dx = d[!(d$SAMPLE_ID %in% t$tip.label), ] %>% filter(RECENT_RUN != TRUE)
# dx2 = dx %>% select(SAMPLE_ID, CURRENT_TAXON, CYTB, SQCL, ddRAD, RECENT_RUN)
# View(dx2)
# #######################
# # all missing inds had too much missing data
# # or couldn't be assigned to CURRENT_TAXON
# #######################