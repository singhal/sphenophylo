rm(list = ls())

library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(ggplot2)
library(cowplot)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
source("scripts/colors.R")

x = read.csv("../Sphenomorphine_speciation/metadata/sphenomorphine_all_v12.csv")
x1 = x %>% group_by(OPERATIONAL_TAXON, CURRENT_TAXON) %>%
  summarize(n = n())
sps = x1$OPERATIONAL_TAXON
names(sps) = x1$CURRENT_TAXON

d = read.csv("data/structure/inds_to_clusters.29June24.csv") %>% 
  select(SAMPLE_ID = ind, pop)
d$species = gsub("_V\\d+$", "", d$pop)
d$OPERATIONAL_TAXON = sps[ d$species ]

# only want to plot more than one sp
d1 = d %>% group_by(OPERATIONAL_TAXON) %>% 
  summarize(n = length(unique(pop))) %>%
  filter(n > 1)
sps = unique(d1$OPERATIONAL_TAXON)

world = ne_countries(scale = "medium", returnclass = "sf")
res = vector("list", length(sps))
names(res) = sps

for (sp in sps) {
  inds = x[which(x$OPERATIONAL_TAXON == sp), ] %>% left_join(d) %>%
    filter(complete.cases(pop))
  res[[sp]] = ggplot() +
    geom_sf(data = world, fill = "white", color = "gray70", size = 0.15) +
    geom_point(data = inds, aes(x = LON, y = LAT, color = pop), size = 2) +
    xlim(112, 154) + ylim(-39, -9) +
    theme_void() + theme(legend.position = "none",
                         plot.title = element_text(face = "italic")) +
    labs(title = gsub("_", " ", sub("[a-z]+_", ". ", sp))) +
    guides(fill = "none", size = "none", alpha = "none", shape = "none") +
    scale_color_manual(values = allcols) +  
    theme(panel.background = element_rect(fill = "aliceblue"),
          panel.border = element_rect(fill = NA, color = "gray30", size = 0.5),
          plot.margin = margin(r = 2, l = 2, t = 0, b = 0))
}

a = plot_grid(plotlist = res[1:25], nrow = 5)
save_plot("figures/cluster_maps1.png", a, base_height = 10, base_width = 10)
b = plot_grid(plotlist = res[26:50], nrow = 5)
save_plot("figures/cluster_maps2.png", b, base_height = 10, base_width = 10)
c = plot_grid(plotlist = res[51:75], nrow = 5)
save_plot("figures/cluster_maps3.png", c, base_height = 10, base_width = 10)
d = plot_grid(plotlist = res[76:86], nrow = 2, ncol = 5)
save_plot("figures/cluster_maps4.png", d, base_height = 4, base_width = 10)
