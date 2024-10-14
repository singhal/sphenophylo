rm(list = ls())
library(dplyr)
library(ape)
library(RColorBrewer)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Ctenotus_speciation/")

d = read.csv("metadata/sphenomorphine_all_v4.csv")
allcols = rep(brewer.pal(12, "Set3"), 30)

# check SqCL tree
t = read.tree("phylogeny/SqCL/concat/concat_ind0.05_loci0.3_all_n5277.phy.treefile.tre")
outs = d[which(d$INGROUP == FALSE), "SAMPLE_ID"]
outs = c(outs, "spheno_unk1", "spheno_unk3")
t = root(t, outs[1], resolve.root = T)
t1 = drop.tip(t, outs)
drop = c("UMMZ_244315_ct_quat", "CUMV_14452_Le_bipe")
t1 = drop.tip(t1, drop)

t1$tip.label[which("CCM5506" == t1$tip.label)] = "NA_CCM5506_Ct_deca"
t1$tip.label[which("CCM6337" == t1$tip.label)] = "NA_CCM6337_Ct_deca"

sps = d[match(t1$tip.label, d$SAMPLE_ID), "CLUSTER_SUBGENUS"]
tips = paste(sps,
             d[match(t1$tip.label, d$SAMPLE_ID), "COLLECTOR_SP"],
             t1$tip.label, sep = " / ")

sps2 = unique(sps)
colors = rep(NA, length(sps))
for (j in 1:length(sps2)) {
  colors[ which(sps == sps2[j]) ] = allcols[j]
}

pdf("~/Desktop/SqCL.pdf", height = 40, width = 10)
par(mar = c(0, 0, 0, 20), xpd = T)
plot(t1, show.tip.label = F)
tiplabels(tips, bg = colors, frame = "rect", adj = c(0, 0.5))
dev.off()

# check ddRAD tree
# Ctenotus
tc = read.tree("~/Dropbox (Personal)/macroevolution/eco_IBD_oz/data/nuc_tree/Ctenotus_nuclear_percent0.6.RAxML.rooted.tre")
sps = d[match(tc$tip.label, d$SAMPLE_ID), "CLUSTER_SUBGENUS"]
tips = paste(sps,
             d[match(tc$tip.label, d$SAMPLE_ID), "COLLECTOR_SP"],
             tc$tip.label, sep = " / ")
sps2 = unique(sps)
colors = rep(NA, length(sps))
for (j in 1:length(sps2)) {
  colors[ which(sps == sps2[j]) ] = allcols[j]
}

pdf("~/Desktop/ddRAD_Ctenotus.pdf", height = 120, width = 10)
par(mar = c(0, 0, 0, 20), xpd = T)
plot(tc, show.tip.label = F)
tiplabels(tips, bg = colors, frame = "rect", adj = c(0, 0.5))
dev.off()

# Lerista
tc = read.tree("~/Dropbox (Personal)/macroevolution/eco_IBD_oz/data/nuc_tree/Lerista_nuclear_percent0.6.RAxML.rooted.tre")
sps = d[match(tc$tip.label, d$SAMPLE_ID), "CLUSTER_SUBGENUS"]
tips = paste(sps,
             d[match(tc$tip.label, d$SAMPLE_ID), "COLLECTOR_SP"],
             tc$tip.label, sep = " / ")
sps2 = unique(sps)
colors = rep(NA, length(sps))
for (j in 1:length(sps2)) {
  colors[ which(sps == sps2[j]) ] = allcols[j]
}

pdf("~/Desktop/ddRAD_Lerista.pdf", height = 80, width = 10)
par(mar = c(0, 0, 0, 20), xpd = T)
plot(tc, show.tip.label = F)
tiplabels(tips, bg = colors, frame = "rect", adj = c(0, 0.5))
dev.off()

# check mtDNA tree
# Ctenotus
tc = read.tree("phylogeny/mtDNA/Ctenotus.treefile")
outs = tc$tip.label[ grep("_Ct_", tc$tip.label, ignore.case = T, invert = T) ]
tc = root(tc, outs[1], resolve.root = T)
tc = drop.tip(tc, outs)
sps = d[match(tc$tip.label, d$SAMPLE_ID), "CLUSTER_SUBGENUS"]
tips = paste(sps,
             d[match(tc$tip.label, d$SAMPLE_ID), "COLLECTOR_SP"],
             tc$tip.label, sep = " / ")
sps2 = unique(sps)
colors = rep(NA, length(sps))
for (j in 1:length(sps2)) {
  colors[ which(sps == sps2[j]) ] = allcols[j]
}
tips2 = paste0(tips, ifelse(d[match(tc$tip.label, d$SAMPLE_ID), "INTROGRESSED"] == "TRUE", "*", ""))
pdf("~/Desktop/mtDNA_Ctenotus.pdf", height = 260, width = 10)
par(mar = c(0, 0, 0, 20), xpd = T)
plot(tc, show.tip.label = F)
tiplabels(tips2, bg = colors, frame = "rect", adj = c(0, 0.5))
dev.off()

# Lerista
tc = read.tree("phylogeny/mtDNA/Lerista.treefile")
outs = tc$tip.label[ grep("_Le_", tc$tip.label, ignore.case = T, invert = T) ]
tc = root(tc, outs[1], resolve.root = T)
tc = drop.tip(tc, outs)
sps = d[match(tc$tip.label, d$SAMPLE_ID), "CLUSTER_SUBGENUS"]
tips = paste(sps,
             d[match(tc$tip.label, d$SAMPLE_ID), "COLLECTOR_SP"],
             tc$tip.label, sep = " / ")
sps2 = unique(sps)
colors = rep(NA, length(sps))
for (j in 1:length(sps2)) {
  colors[ which(sps == sps2[j]) ] = allcols[j]
}
tips2 = paste0(tips, ifelse(d[match(tc$tip.label, d$SAMPLE_ID), "INTROGRESSED"] == "TRUE", "*", ""))
pdf("~/Desktop/mtDNA_Lerista.pdf", height = 80, width = 10)
par(mar = c(0, 0, 0, 20), xpd = T)
plot(tc, show.tip.label = F)
tiplabels(tips2, bg = colors, frame = "rect", adj = c(0, 0.5))
dev.off()
