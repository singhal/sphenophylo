rm(list = ls())
library(tidyverse)
library(cowplot)
library(ape)

rearrange_species_names <- function(df, scol1, scol2) {
  species1 = mapply(function(x, y) min(x, y), df[, scol1], df[, scol2])
  species2 = mapply(function(x, y) max(x, y), df[, scol1], df[, scol2])
  df[,scol1] = species1
  df[,scol2] = species2
  return(df)
}

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")

source("scripts/colors.R")
names(genera) = NULL

d = read.csv("../Sphenomorphine_speciation/metadata/sphenomorphine_all_v13.csv")
x = read.csv("data/species_to_OTU/mapping.amended.csv") %>% filter(complete.cases(delimitation))
t = read.tree("../Sphenomorphine_speciation/phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")


# get Fst, div time values for species
source("../Sphenomorphine_speciation/scripts/create_master_dataset.R")
# keeping all sisters
# ((A, B), C) - C is sister to A and B; A is only sister to B
dd1 = dd %>% filter(sister_species == TRUE)

# fix naming to operational
d1 = d %>% group_by(CURRENT_TAXON, OPERATIONAL_TAXON) %>% summarize(n = n()) %>%
  dplyr::select(CURRENT_TAXON, OPERATIONAL_TAXON)
namefix = d1$OPERATIONAL_TAXON
names(namefix) = d1$CURRENT_TAXON

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
y = read.csv("data/structure/inds_to_clusters.29June24.csv")
y$species = gsub("_V\\d+", "", y$pop)
y$OPERATIONAL_TAXON = namefix[ y$species ]

dd1$species1 = namefix[ dd1$species1 ]
dd1$species2 = namefix[ dd1$species2 ]
dd1$type1 = x[match(dd1$species1, x$OPERATIONAL_TAXON), "delimitation"]
dd1$type2 = x[match(dd1$species2, x$OPERATIONAL_TAXON), "delimitation"]

# okay to not do anything with combined?
dd1$type = ifelse(dd1$type1 == "split" & dd1$type2 == "split", "split", "other")
# relationships with musivus are undefined
dd1[dd1$species2 == "Eremiascincus_musivus", "type"] = "other"

###################
# first make a figure to show 
# why we focus on just a few
# types of comparisons
###################

t1 = read.tree("data/diversification/species_level_phylogeny.OPERATIONAL_TAXON.tre")
t1 = cophenetic.phylo(t1) %>% as.data.frame() %>% 
  rownames_to_column("species1") %>%
  gather(species2, phydist, -species1)

y1 = dd1 %>% select(species1, species2, fst) %>% 
  mutate(taxonomy = "operational") %>% left_join(t1)

t2 = read.tree("data/diversification/species_level_phylogeny.MORPHOLOGICAL_TAXON.tre")
t2 = cophenetic.phylo(t2) %>% as.data.frame() %>% 
  rownames_to_column("species1") %>%
  gather(species2, phydist, -species1)

y2 = list.files("data/divergence_morph/", full.names = T)
y2 = do.call("rbind", lapply(y2, read.csv)) %>%
  filter(metric == "fst", taxon_type == "fst") %>% select(species1, species2, fst = submetric)
y2 = y2 %>% mutate(taxonomy = "morphological") %>% left_join(t2)

t3 = read.tree("data/diversification/species_level_phylogeny.THRESHOLD_TAXON.tre")
t3 = cophenetic.phylo(t3) %>% as.data.frame() %>% 
  rownames_to_column("species1") %>%
  gather(species2, phydist, -species1)

y3 = list.files("data/divergence_threshold/", full.names = T)
y3 = do.call("rbind", lapply(y3, read.csv)) %>%
  filter(metric == "fst", taxon_type == "fst") %>% select(species1, species2, fst = submetric)
y3 = y3 %>% mutate(taxonomy = "threshold") %>% left_join(t3)

t4 = read.tree("data/diversification/species_level_phylogeny.CLUSTERS.tre")
t4 = cophenetic.phylo(t4) %>% as.data.frame() %>% 
  rownames_to_column("species1") %>%
  gather(species2, phydist, -species1)

y4 = list.files("data/divergence_incipient/structure//", full.names = T, pattern = ".csv")
y4 = do.call("rbind", lapply(y4, read.csv)) %>%
  filter(metric == "fst", submetric == "fst") %>% select(species1, species2, fst = value)
y4 = y4 %>% mutate(taxonomy = "incipient") %>% left_join(t4)

all = rbind(y1, y2, y3, y4)
all$btime = all$phydist / 2

a = ggplot(all) +
  geom_density(aes(fst, fill = taxonomy), alpha = 0.3, linewidth = 0.1) +
  scale_fill_manual(values = tcols) +
  xlab(expression(F[ST])) +
  theme_classic()
b = ggplot(all) +
  geom_density(aes(btime, fill = taxonomy), alpha = 0.3, linewidth = 0.1) +
  scale_fill_manual(values = tcols) +
  xlab("branching time") +
  theme_classic()

prow <- plot_grid(
  a + theme(legend.position="none"),
  b + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)

legend <- get_legend(
  # create some space to the left of the legend
  a + theme(legend.box.margin = margin(0, 0, 0, 0))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
alldiv = plot_grid(prow, legend, rel_widths = c(3, .4))
save_plot("figures/all_comparisons.png", alldiv, base_height = 3, base_width = 9)

###################
# LTT 
###################

outs = d[which(d$INGROUP == FALSE), "SAMPLE_ID"]
t1 = drop.tip(t, outs)
t1$node.label = seq(Ntip(t1) + 1, Ntip(t1) + Nnode(t1))
times = branching.times(t1)

tt = data.frame(
  node = t1$node.label,
  age = times[as.character(t1$node.label)]
)
tt[rev(order(tt$age)),  "ntip"] = seq(2, Ntip(t1)) 
maxht = max(phytools::nodeHeights(t1))
tt$time = -1 * tt$age

tt$type = "all"
# # identify within
# for (i in 1:nrow(tt)) {
#   node = tt$node[i]
#   desc = phytools::getDescendants(t1, node)
#   desc = t1$tip.label[ desc[desc < (Ntip(t1) + 1)] ]
#   sps = unique(d[d$SAMPLE_ID %in% desc, "CURRENT_TAXON"])
#   if (length(sps) == 1) {
#     tt[tt$node == node, "type"] = "within"
#   }
# }
# sps = unique(d$CURRENT_TAXON)
# sps = sps[complete.cases(sps)]
# # identify crown ages & stem (sister) ages
# for (sp in sps) {
#   desc = d[which(d$CURRENT_TAXON == sp), "SAMPLE_ID"]
#   if (length(desc) == 1) {
#     parent = phytools::getParent(t1, which(t$tip.label == desc))
#     tt[tt$node == parent, "type"] = "stem"
#   } else {
#     node = phytools::findMRCA(t1, desc)
#     parent = phytools::getParent(t1, node)
#     # cat(sp, "\n")
#     tt[tt$node == node, "type"] = "crown"
#     tt[tt$node == parent, "type"] = "stem"
#   }
# }
# 
# table(tt$type)

pops = unique(y$pop)
for (i in 1:length(pops)) {
  inds = y[y$pop == pops[i], "ind"]
  keynode = phytools::findMRCA(t1, inds)
  tt[tt$node == keynode, "type"] = "incipient"
}

aa = ggplot(tt, aes(time, ntip)) +
  geom_point(col = "gray80")  +
  geom_point(data = tt %>% filter(type == "incipient"), 
             aes(time, ntip),
             fill = genera[1], pch = 21)  +
  geom_point(data = tt %>% filter(type == "incipient"), 
             aes(time, ntip),
             col = genera[1], pch = 16)  +
  xlab("time (myr)") +
  ylab("# of tips") +
  theme_classic()

###################
# divtime 
###################
combos = function(x) {
  npop = unique(x$pop)
  if (length(npop) == 1) {
    matrix(NA, nrow = 1, ncol = 2)
  } else if (length(npop) == 2) {
    return(matrix(npop, nrow = 1, ncol = 2))
  } else {
    return(t(combn(npop, 2)))
  }
}

y1 = split(y, y$species)
y2 = as.data.frame(do.call("rbind", lapply(y1, combos)))
names(y2) = c("species1", "species2")
y3 = y2 %>% filter(complete.cases(species1)) %>% 
  mutate(type = "incipient", divtime = NA)

tt = branching.times(t)
for (i in 1:nrow(y3)) {
 sps = y3[i, c("species1", "species2")] 
 inds =  y[y$pop %in% sps, "ind"]
 inds = inds[inds %in% t$tip.label]
 if (length(inds) > 1) {
   y3$divtime[i] = tt[ phytools::findMRCA(t, inds) - Ntip(t) ]
 }
}
y4 = rearrange_species_names(y3, "species1", "species2")

###################
# Fst 
###################

# get Fst for phylogroups
f = list.files("data/divergence_incipient/structure/", full.names = T)
f1 = do.call("rbind", lapply(f, read.csv)) %>% filter(metric == "fst")
f2 = f1 %>% spread(submetric, value) %>%
  filter(fst_counts > 1000) %>% 
  dplyr::select(-metric)
f3 = rearrange_species_names(f2, "species1", "species2")

# ###################
# # migration etc 
# ###################
 
m = read.csv("data/demography/structure/demography_summary.18June24.csv")
m1 = rearrange_species_names(m, "species1", "species2") %>% 
   dplyr::select(species1, species2, M12, M21, T)
 
y5 = left_join(y4, f3) %>% left_join(m1)
dd2 = rbind(dd1 %>% dplyr::select(species1, species2, divtime = phydist / 2, fst, type, M12, M21, T),
             y5 %>% dplyr::select(species1, species2, divtime, fst, type, M12, M21, T)) %>%
   mutate(rank_fst = dense_rank(fst)) %>% rowwise() %>%
   mutate(xM = max(M12, M21)) %>% ungroup() %>% mutate(rank_M = dense_rank(xM))

###################
# plots
###################
tcols2 = tcols
names(tcols2) = c('incipient', 'split', "other")

a = ggplot(dd2) + 
  geom_density(aes(x = divtime, fill = type), alpha = 0.3, linewidth = 0.1) + 
  scale_fill_manual(values = tcols2) +
  xlab("divergence time (Mya)") +
  theme_classic() +
  theme(legend.position = c(0.8, 0.8)) +
  theme(legend.title=element_blank())
  
# ggplot(dd2) +
#   geom_density(aes(x = fst, fill = type), alpha = 0.3, linewidth = 0.1) + 
#   scale_fill_manual(values = genera) +
#   xlab(expression(F[ST])) +
#   theme_classic()
# ggplot(dd2) + 
#   geom_density(aes(x = xM, fill = type), alpha = 0.3, linewidth = 0.1) + 
#   scale_fill_manual(values = genera) +
#   xlab("max migration rate (M)") +
#   scale_x_log10() +
#   theme_classic()

b = ggplot(dd2) +
  geom_point(aes(x = rank_fst, y = fst), color = "gray", shape = 16) +
  # geom_point(data = dd2 %>% filter(type == "split"), aes(x = rank_fst, y = fst), fill = tcols2["split"], shape = 21) +
  # geom_point(data = dd2 %>% filter(type == "split"), aes(x = rank_fst, y = fst), color = tcols2["split"], shape = 16) +
  geom_point(data = dd2 %>% filter(type == "incipient"), aes(x = rank_fst, y = fst), fill = tcols2["incipient"], shape = 21) +
  geom_point(data = dd2 %>% filter(type == "incipient"), aes(x = rank_fst, y = fst), color = tcols2["incipient"], shape = 16) +
  xlab(expression("species comparisons (ranked by" ~ F[ST] ~ ")")) +
  ylab(expression(F[ST])) +
  theme_classic() +
  theme(axis.ticks.x=element_blank())
c = ggplot(dd2, aes(fst, xM +
                      0.001, col = type)) + 
  geom_point() + theme_classic() + scale_y_log10() +
  scale_color_manual(values = c(genera[1], genera[3], genera[2])) + 
  xlab(expression(F[ST])) + ylab("max migration rate (M)") +
  theme(legend.title = element_blank())
# ggplot(dd2) +
#   geom_point(aes(x = rank_M, y = xM), color = "gray", shape = 16) +
#   geom_point(data = dd2 %>% filter(type == "within-species"), aes(x = rank_M, y = xM), fill = "green", shape = 21) +
#   geom_point(data = dd2 %>% filter(type == "within-species"), aes(x = rank_M, y = xM), color = "green", shape = 16) +
#   xlab("") +
#   ylab(expression(F[ST])) +
#   theme_classic() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())

prow <- plot_grid(
  aa,
  a,
  b + theme(legend.position="none"),
  # c + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)

save_plot("figures/clusters_species.png", prow, 
          base_height = 3, base_width = 9)

# fit curve
log_fit <- nls(fst ~ a + b * log(divtime), data = dd2, start = list(a = 0, b = 1))
# Generate points for the fitted curve
curve_data <- data.frame(divtime = seq(min(dd2$divtime), max(dd2$divtime), length.out = 100))
curve_data$fst = predict(log_fit, newdata = curve_data)

fstplt = ggplot(dd2, aes(divtime, fst)) +
  geom_point(shape = 16, alpha = 0.5, aes(color = type)) +
  geom_line(data = curve_data, aes(x = divtime, y = fst), color = "red", linetype = "dashed") +
  # geom_vline(xintercept = 2, col = "red") +
  theme_classic() +
  scale_color_manual(values = tcols2) + 
  xlab("divergence time") + 
  ylab(expression(F[ST]))
save_plot("figures/fst_divtime.png", fstplt, base_height = 3, base_width = 4)
