rm(list = ls())
library(dplyr)
library(ape)
library(RColorBrewer)
library(reshape2)

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

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
d = read.csv("../Sphenomorphine_speciation/metadata/sphenomorphine_all_v13.csv")

# these individuals have ddRAD data but very low quality
dropinds = c("SAMR_56616_Ct_atla", "WAMR_135958_Ct_pant",
             "WAMR_98141_Le_pete", "SAMR_57650_Ct_atla")
# all outgroups are monophyletic
# except for C. burb -- drop
dropinds = c(dropinds, d[which(d$CURRENT_TAXON == "Ctenotus_burbidgei"), "SAMPLE_ID"])

sqcl = read.tree("../Sphenomorphine_speciation/phylogeny/SqCL/concat/concat_ind0.05_loci0.3_all_n5277.rooted.tre")
# restree = read.tree("../Sphenomorphine_speciation/phylogeny/dating/synthetic_tree.all.rooted.18April2023.tre")

pp = list.files('../Sphenomorphine_speciation/phylogeny/synthetic_trees/',
     pattern = "all.constraint.treefile",
     full.names = T)
res = vector("list", length(pp))
slopes = rep(NA, length(pp))
for (i in 1:length(pp)) {
    x = read.tree(pp[i])

    # get rid of bad quality individuals
    x = drop.tip(x, dropinds)

    x$tip.label = gsub(".outgroup", "", x$tip.label)
    shared = intersect(x$tip.label, sqcl$tip.label)

    tt1 = keep.tip(x, shared)
    tt2 = keep.tip(sqcl, shared)

    pd = get_phy_dist(list(tt1, tt2))
    
    slopes[i] = coef(lm(pd$ddrad ~ pd$sqcl))[2]
    pltix = gsub(".all.constraint.*", "", pp[i])
    pltix = gsub(".*//", "", pltix)
    pd$group = pltix
    pd$nind = length(shared)
    res[[i]] = pd
}

res2 = do.call("rbind", res)

a = ggplot(res2 %>% filter(nind > 3), aes(sqcl, ddrad)) + 
  geom_point(alpha = 0.7, shape = 16) +
  geom_smooth(method='lm', col = "red", se=F) +
  xlab("target capture tree distance") + ylab("synthetic tree distance") +
  facet_wrap(~group, nrow = 8) + theme_classic() + 
  theme(panel.border = element_rect(fill = NA, color = "gray30", size = 0.5),
        strip.text = element_text(size = 9, face = "italic", 
                                  margin = margin(t = 0, r = 0, b = 2, l = 0)), 
        strip.background = element_blank())
save_plot("figures/age_correlations.png", a, base_height = 12, 
          base_width = 8)

