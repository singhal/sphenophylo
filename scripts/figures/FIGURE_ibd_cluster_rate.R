rm(list = ls())
library(ape)
library(sf)
library(tidyverse)
library(nlme)
library(vegan)
library(BAMMtools)

source("scripts/colors.R")

world_map <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
t = read.tree("../Sphenomorphine_speciation/phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")
dd = cophenetic.phylo(t)

d = read.csv("../Sphenomorphine_speciation/metadata/sphenomorphine_all_v11.csv")
d = d[complete.cases(d$CURRENT_TAXON), ]
sps = table(d$CURRENT_TAXON)
# these are the units that have enough sampling
sps = names(sps[sps > 3])

#########################
# calculate IBD
#########################
ibdres = vector("list", length(sps))
names(ibdres) = sps
for (sp in sps) {
  # sp = "Ctenotus_pantherinus"
  dx = d %>% filter(CURRENT_TAXON == sp) %>% filter(complete.cases(LAT))
  inds = dx$SAMPLE_ID
  
  # phylogenetic distances
  phydist = dd[inds, inds]

  # geographic distances
  pts = dx %>% st_as_sf(coords = c("LON", "LAT")) %>% 
    st_set_crs(sf::st_crs(world_map))
  geodist = matrix(as.numeric(st_distance(pts)), nrow = length(inds), ncol = length(inds))
  dimnames(geodist) = list(inds, inds)
  loggeodist = log(geodist + 1)
  
  # do mantel test
  mantelr = mantel(phydist, loggeodist)
  
  # get slope
  xy <- t(combn(inds, 2))
  df = data.frame(xy, loggeodist=loggeodist[xy], phydist = phydist[xy])
  dlm = lm(df$phydist ~ df$loggeodist)
  
  ibdres[[sp]] = c(sp, length(inds), coef(dlm)[2], mantelr$signif, mantelr$statistic)
}

res = as.data.frame(do.call("rbind", ibdres))
names(res) = c("species", "nInds", "beta_ibd", "signif", "r")
for (i in 2:5) {
  res[,i] = as.numeric(res[,i])
}

#########################
# get # of clusters
#########################

res$cluster = NA
for (sp in sps) {
  file = paste0("data/structure/fastStructure/", sp, ".chooseK.csv")
  if (file.exists(file)) {
  c = read.csv(file)
  res[res$species == sp, "cluster"] = c[2, "K"]
  }
}

#########################
# get spec rates
#########################

tc = read.tree("data/diversification/species_level_phylogeny.CURRENT_TAXON.tre")
load("data/diversification/species_level_phylogeny.CURRENT_TAXON.CLaDS.Rdata")
cc = CladsOutput
rates = cc$lambdatip_map
names(rates) = tc$tip.label
res$lambda = rates[res$species]

res$genus = gsub("_\\S+", "", res$species)
res$genus = ifelse(res$genus %in% c("Ctenotus", "Lerista"), res$genus, "other")
res$log_ibd = log(res$beta_ibd)

#######################
# plots
#######################

res3 = res %>% filter(complete.cases(cluster)) %>% group_by(cluster) %>%
  summarize(beta_ibd = mean(beta_ibd))
a = ggplot(res %>% filter(complete.cases(cluster)), aes(beta_ibd, cluster)) +
  geom_jitter(shape = 21, aes(fill = genus), height = 0.1, alpha = 0.7) +
  geom_crossbar(data=res3, aes(xmin = beta_ibd, xmax = beta_ibd),
                size=0.5,col="black", width = .5) +
  theme_classic() + xlab(expression(beta[IBD])) +
  ylab("# of clusters") + scale_x_log10() +
  scale_fill_manual(values = genera)
b = ggplot(res, aes(beta_ibd, lambda)) +
  geom_point(shape = 21, aes(fill = genus)) +
  theme_classic() + xlab(expression(beta[IBD])) +
  ylab(expression(lambda[CLaDS])) + scale_x_log10() +
  scale_fill_manual(values = genera)
res2 = res %>% filter(complete.cases(cluster)) %>% group_by(cluster) %>%
  summarize(lambda = mean(lambda))
c = ggplot(res %>% filter(complete.cases(cluster)), aes(as.factor(cluster), lambda)) +
  geom_jitter(width = 0.1, pch = 21, aes(fill= genus), alpha = 0.7) +
  geom_crossbar(data=res2, aes(ymin = lambda, ymax = lambda),
                size=0.5,col="black", width = .5) +
  theme_classic() + xlab("# of clusters") +
  ylab(expression(lambda[CLaDS])) +
  scale_fill_manual(values = genera)

prow <- cowplot::plot_grid(
  a + theme(legend.position="none"),
  b + theme(legend.position="none"),
  c + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)
prow

ll <- cowplot::get_legend(
  # create some space to the left of the legend
  b + theme(legend.box.margin = margin(0, 0, 0, 0))
)

# add the legend to the row we made earlier
corr = cowplot::plot_grid(prow, ll, rel_widths = c(3, 0.3))
cowplot::save_plot("figures/correlations.png", corr, base_width = 9, 
                   base_height = 3)

####################
# phylo linear regressions
####################
# ibd to lambda
res2 = res %>% filter(complete.cases(log_ibd))
rownames(res2) = res2$species
t2 = keep.tip(tc, res2$species)
res3 = res2[t2$tip.label, ]
x = gls(lambda ~ log_ibd, 
           correlation = corBrownian(phy = t2),
           data = res3, method = "ML")
summary(x)

# ibd to cluster
res2 = res %>% filter(complete.cases(log_ibd), complete.cases(cluster))
rownames(res2) = res2$species
t2 = keep.tip(tc, res2$species)
res3 = res2[t2$tip.label, ]
y = gls(lambda ~ cluster, 
        correlation = corBrownian(phy = t2),
        data = res3, method = "ML")
summary(y)

# cluster, speciation
res2 = res %>% filter(complete.cases(cluster))
rownames(res2) = res2$species
t2 = keep.tip(tc, res2$species)
res3 = res2[t2$tip.label, ]
z = gls(lambda ~ cluster, 
        correlation = corBrownian(phy = t2),
        data = res3, method = "ML")
summary(z)

