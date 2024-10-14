rm(list = ls())
library(tidyverse)
library(cowplot)
library(ape)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(ggtree)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")

source("scripts/colors.R")

world = ne_countries(scale = "medium", returnclass = "sf")
d = read.csv("../Sphenomorphine_speciation/metadata/sphenomorphine_all_v13.csv")
dx = d
t = read.tree("../Sphenomorphine_speciation/phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")
x = read.csv("data/structure/inds_to_clusters.29June24.csv")
y = read.csv("data/threshold_taxonomy/threshold_taxa.csv")

lins = c("Ctenotus_schomburgkii_1", "Ctenotus_schomburgkii_2", "Ctenotus_schomburgkii_3", "Ctenotus_kutjupa")
inds = d[d$CURRENT_TAXON %in% lins, "SAMPLE_ID"]

tx = extract.clade(t, findMRCA(t, inds))
dx = d[d$SAMPLE_ID %in% inds, ]
dy = dx[which(dx$ddRAD == TRUE), ]
dy = dy[which(dy$RECENT_RUN == FALSE), ]
dy$CLUSTER = coalesce(x[match(dy$SAMPLE_ID, x$ind), "pop"], dy$CURRENT_TAXON)
ty = keep.tip(tx, dy$SAMPLE_ID)

###################
# do the morphological
###################

a = ggtree(tx, ladderize = F, color = allcols[1], linewidth = 0.3) +
  geom_cladelab(node = Ntip(tx) + 1, label = "C. schomburgkii", 
                align = TRUE, offset = 0.2, 
                textcolor = allcols[1], barcolor = allcols[1]) + 
  xlim(0, 15) 

b = ggplot() +
  geom_sf(data = world, fill = "white", color = "gray70", size = 0.15) +
  geom_point(data = dx, aes(x = LON, y = LAT), size = 2, fill = allcols[1], shape = 21) +
  xlim(112, 154) + ylim(-39, -9) +
  theme_void() + theme(legend.position = "none",
                       plot.title = element_text(face = "italic")) +
  guides(fill = "none", size = "none", alpha = "none", shape = "none") +
  theme(panel.background = element_rect(fill = "aliceblue"),
        panel.border = element_rect(fill = NA, color = "gray30", size = 0.5),
        plot.margin = margin(r = 2, l = 2, t = 0, b = 0))


###################
# do the operational
###################

lincols = allcols[1:4]
names(lincols) = lins

# color branches by clade
gps = lapply(lins, function(sp) dx[which(dx$CURRENT_TAXON == sp), "SAMPLE_ID"])
names(gps) = lins
tx2 = groupOTU(tx, gps)

cols1 = c("black", lincols)
names(cols1)[1] = "0"

c = ggtree(tx2, linewidth = 0.3, ladderize = FALSE, aes(color = group)) +
  scale_color_manual(values = cols1) +
  theme(legend.position = "none")
for (i in 1:length(lins)) {
  lintips = dx[dx$CURRENT_TAXON == lins[i], "SAMPLE_ID"]
  nodex = findMRCA(tx, tips = lintips)
  linname = gsub("_", " ", gsub("Ctenotus_", "C. ", lins[i]))
  linname = gsub("schomburgkii", "scho.", linname)
  c = c + geom_cladelab(node = nodex, label = linname, 
                align = TRUE, offset = 0.2, 
                textcolor = lincols[i], barcolor = lincols[i])
}
c = c + xlim(0, 15) 

d = ggplot() +
  geom_sf(data = world, fill = "white", color = "gray70", size = 0.15) +
  geom_point(data = dx, aes(x = LON, y = LAT, fill = CURRENT_TAXON), size = 2, shape = 21) +
  xlim(112, 154) + ylim(-39, -9) +
  theme_void() + theme(legend.position = "none") +
  scale_fill_manual(values = lincols) +
  guides(fill = "none", size = "none", alpha = "none", shape = "none") +
  theme(panel.background = element_rect(fill = "aliceblue"),
        panel.border = element_rect(fill = NA, color = "gray30", size = 0.5),
        plot.margin = margin(r = 2, l = 2, t = 0, b = 0))


###################
# do the threshold
###################

y1 = y[grep("scho", y$CURRENT_TAXON), ]
y2 = y[grep("kutj", y$CURRENT_TAXON), ]
y1 = rbind(y1, y2)
lins = unique(y1$sp)

lincols = allcols[1:length(lins)]
names(lincols) = lins

# color branches by clade
gps = lapply(lins, function(sp) y1[which(y1$sp == sp), "ind"])
names(gps) = lins
tx2 = groupOTU(tx, gps)

cols1 = c("black", lincols)
names(cols1)[1] = "0"

g = ggtree(tx2, linewidth = 0.3, ladderize = FALSE, aes(color = group)) +
  scale_color_manual(values = cols1) +
  theme(legend.position = "none")
for (i in 1:length(lins)) {
  lintips = y1[which(y1$sp == lins[i]), "ind"]
  nodex = findMRCA(tx, tips = lintips)
  linname = gsub("sp", "C. sp", lins[i])
  g = g + geom_cladelab(node = nodex, label = linname, 
                        align = TRUE, offset = 0.2, 
                        textcolor = lincols[i], barcolor = lincols[i])
}
g = g + xlim(0, 15) 

y1 = left_join(y1, dx, by = c("ind" = "SAMPLE_ID"))

h = ggplot() +
  geom_sf(data = world, fill = "white", color = "gray70", size = 0.15) +
  geom_point(data = y1, aes(x = LON, y = LAT, fill = sp), size = 2, shape = 21) +
  xlim(112, 154) + ylim(-39, -9) +
  theme_void() + theme(legend.position = "none") +
  scale_fill_manual(values = lincols) +
  guides(fill = "none", size = "none", alpha = "none", shape = "none") +
  theme(panel.background = element_rect(fill = "aliceblue"),
        panel.border = element_rect(fill = NA, color = "gray30", size = 0.5),
        plot.margin = margin(r = 2, l = 2, t = 0, b = 0))


###################
# do the incipient
###################

cls = unique(dy$CLUSTER)
clcols = allcols[1:length(cls)]
names(clcols) = cls

# color branches by clade
gps = lapply(cls, function(sp) dy[which(dy$CLUSTER == sp), "SAMPLE_ID"])
names(gps) = cls
ty2 = groupOTU(ty, gps)

cols1 = c("black", clcols)
names(cols1)[1] = "0"

e = ggtree(ty2, linewidth = 0.3, ladderize = FALSE, aes(color = group)) +
  scale_color_manual(values = cols1) +
  theme(legend.position = "none")
for (i in 1:length(cls)) {
  lintips = dy[dy$CLUSTER == cls[i], "SAMPLE_ID"]
  if (length(lintips) > 1) {
    nodex = findMRCA(ty, tips = lintips)
  } else {
    nodex = which(ty$tip.label == lintips)
  }
  linname = gsub("_", " ", gsub("Ctenotus_", "C. ", cls[i]))
  linname = gsub("schomburgkii", "scho.", linname)
  linname = gsub("V", "lin. ", linname)
  e = e + geom_cladelab(node = nodex, label = linname, 
                        align = TRUE, offset = 0.2, 
                        textcolor = clcols[i], barcolor = clcols[i])
}
e = e + xlim(0, 15) 

f = ggplot() +
  geom_sf(data = world, fill = "white", color = "gray70", size = 0.15) +
  geom_point(data = dy, aes(x = LON, y = LAT, fill = CLUSTER), size = 2, shape = 21) +
  xlim(112, 154) + ylim(-39, -9) +
  theme_void() + theme(legend.position = "none") +
  scale_fill_manual(values = clcols) +
  guides(fill = "none", size = "none", alpha = "none", shape = "none") +
  theme(panel.background = element_rect(fill = "aliceblue"),
        panel.border = element_rect(fill = NA, color = "gray30", size = 0.5),
        plot.margin = margin(r = 2, l = 2, t = 0, b = 0))

xx = a + c + e + g + b + d + f + h + plot_layout(ncol = 4)
cowplot::save_plot("figures/explanatory_figure.pdf", xx, base_height = 6, base_width = 13)
