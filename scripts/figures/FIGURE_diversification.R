rm(list = ls())
library(ggplot2)
library(tidyverse)
library(ape)
library(phytools)
library(BAMMtools)
library(patchwork)
library(ggtree)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
source("scripts/colors.R")
source("scripts/DR_functions.R")


tc = read.tree("data/diversification/species_level_phylogeny.OPERATIONAL_TAXON.tre")
tn = read.tree("data/diversification/species_level_phylogeny.MORPHOLOGICAL_TAXON.tre")
tl = read.tree("data/diversification/species_level_phylogeny.CLUSTERS.tre")
tt = read.tree("data/diversification/species_level_phylogeny.THRESHOLD_TAXON.tre")

#################
# LTT plot
#################

get_ltt <- function(t) {
  t$node.label = seq(Ntip(t) + 1, Ntip(t) + Nnode(t))
  times = branching.times(t)
  
  tt = data.frame(
    node = t$node.label,
    age = times[as.character(t$node.label)]
  )
  tt[rev(order(tt$age)),  "ntip"] = seq(2, Ntip(t)) 
  maxht = max(nodeHeights(t))
  tt$time = -1 * tt$age
  tt = rbind(tt, c(NA, NA, Ntip(t), 0))
  return(tt)
}

ltt1 = get_ltt(tc)
ltt1$type = "operational"
ltt2 = get_ltt(tn)
ltt2$type = "morphological"
ltt3 = get_ltt(tl)
ltt3$type = "incipient"
ltt4 = get_ltt(tt)
ltt4$type = "threshold"
ltt = rbind(ltt1, ltt2, ltt3, ltt4)

a = ggplot(ltt, aes(time, ntip)) +
  geom_line(aes(col = type)) +
  ylab("# of lineages") +
  scale_y_log10() +
  scale_color_manual(values = tcols) +
  theme_classic() +
  annotate("text", y=310, x=-16, label= "morphological", col = tcols["morphological"]) +
  annotate("text", y=220, x=-16, label= "operational", col = tcols["operational"]) +
  annotate("text", y=158, x=-16, label= "incipient", col = tcols["incipient"]) +
  annotate("text", y=115, x=-16, label= "threshold", col = tcols["threshold"]) +
  theme(legend.position = "none")

########################
# CLaDS output AND DR
########################

load("data/diversification/species_level_phylogeny.OPERATIONAL_TAXON.CLaDS.Rdata")
cc = CladsOutput
r1 = data.frame(rate = cc$lambdatip_map,
                tip = tc$tip.label, 
                type = "operational")

load("data/diversification/species_level_phylogeny.MORPHOLOGICAL_TAXON.CLaDS.Rdata")
cn = CladsOutput
r2 = data.frame(rate = cn$lambdatip_map,
                tip = tn$tip.label, 
                type = "morphological")

load("data/diversification/species_level_phylogeny.CLUSTERS.CLaDS.Rdata")
cl = CladsOutput
r3 = data.frame(rate = cl$lambdatip_map,
                tip = tl$tip.label, 
                type = "incipient")

load("data/diversification/species_level_phylogeny.THRESHOLD_TAXON.CLaDS.Rdata")
cl = CladsOutput
r4 = data.frame(rate = cl$lambdatip_map,
                tip = tt$tip.label, 
                type = "threshold")


r5 = rbind(r1, r2, r3)
r5$inference = "CLaDS"

r4 = rbind(data.frame(rate = DR(tc), type = "operational"),
      data.frame(rate = DR(tn), type = "morphological"),
      data.frame(rate = DR(tl), type = "incipient"),
      data.frame(rate = DR(tt), type = "threshold"))
r4$inference = "DR statistic"
r4$tip = rownames(r4)

rates = rbind(r4, r5)
rates$genus = gsub("_\\S+", "", rates$tip)
rates$genus2 = ifelse(rates$genus %in% c("Ctenotus", "Lerista"), rates$genus, "other")

rates$type2 = factor(rates$type, 
                     levels = c("morphological", "operational", "incipient", "threshold"))

# https://wilkelab.org/ggtext/
b = ggplot(rates, aes(genus2, rate, fill = type2)) +
  geom_boxplot(outlier.size = 0.2) + 
  xlab("") + ylab("speciation rate") +
  theme_classic() +
  scale_fill_manual(values = tcols) +
  theme(legend.title = element_blank()) +
  facet_wrap(~inference, scales="free_y") +
  theme(panel.border = element_rect(fill = NA, color = "gray30", linewidth = 0.5),
        strip.text = element_text(size = 9,margin = margin(t = 0, r = 0, b = 2, l = 0)), 
        strip.background = element_blank())

#################
# BAMM output
#################

############ 
# operational
#############

edc = getEventData(tc, 
                    eventdata = "data/diversification/OPERATIONAL_TAXON_event_data.txt",
                    burnin=0.25, nsamples = 1000)
shifts_c =  credibleShiftSet(edc, 
                        expectedNumberOfShifts=1, 
                        threshold=3, set.limit = 0.93)

# plotRateThroughTime(edc, ratetype="speciation")

# get most likely shift
keynode1 = shifts_c$shiftnodes[[1]]
# get probability 
shifts_c$frequency[1]
r1 = getTipRates(edc)
r1 = data.frame(rate = r1$lambda.avg,
                tip = names(r1$lambda.avg), 
                type = "operational")

rtt1 <- getRateThroughTimeMatrix(edc, node = keynode1, 
                                 nodetype = "include")
rtt2 <- getRateThroughTimeMatrix(edc, node = keynode1, 
                                 nodetype = "exclude")

pdf("figures/slowdown.pdf", width = 4, height = 8)
par(mfrow = c(2, 1))
plot.bammdata(edc, lwd=2, legend = T, breaksmethod='linear',
              logcolor = T)
axisPhylo(backward = F)

## Speciation quantiles: plot 90% CIs
qq1 <- apply(rtt1$lambda, 2, quantile, c(0.05, 0.5, 0.95))
xv1 <- c(rtt1$times, rev(rtt1$times))
yv1 <- c(qq1[1,], rev(qq1[3,]))

plot.new()
plot.window(xlim=c(0, max(rtt1$times)), ylim=c(0, max(qq1)))

## Add the confidence polygon on rate distributions:
polygon(xv1, yv1, col=alpha(tcols[3], 0.3), border=FALSE)
## Add the median rate line:
lines(rtt1$times, qq1[2,], lwd=3, col=tcols[3])

## Speciation quantiles: plot 90% CIs
qq2 <- apply(rtt2$lambda, 2, quantile, c(0.05, 0.5, 0.95))
xv2 <- c(rtt2$times, rev(rtt2$times))
yv2 <- c(qq2[1,], rev(qq2[3,]))
## Add the confidence polygon on rate distributions:
polygon(xv2, yv2, col=alpha(tcols[1], 0.3), border=FALSE)
## Add the median rate line:
lines(rtt2$times, qq2[2,], lwd=3, col=tcols[1])

## Add axes
axis(1, at=seq(-5, 35, by=5))
mtext("time (myr)", side = 1, line = 2)
axis(2, at=seq(-0.5, 2, by=0.5), las=1)
mtext("speciation rate", side = 2, line = 2.5)
dev.off()



t1 = ggtree(tc, ladderize = F, color = tcols["operational"], linewidth = 0.3) +
  geom_nodepoint(aes(subset = node %in% keynode1),
                 alpha = 0.5, size = 5)
for (i in c("Ctenotus", "Lerista")) {
  cnode = phytools::findMRCA(tc, tc$tip.label[ grep(i, tc$tip.label) ])
  t1 = t1 + geom_cladelab(node=cnode, label = i, 
                          align=TRUE, 
                          offset = .2,
                          textcolor = tcols["operational"],
                          barcolor =  tcols["operational"]) 
}
t1 = t1 + xlim(0, 30)

############ 
# recognized
#############

edn = getEventData(tn, 
                   eventdata = "data/diversification/MORPHOLOGICAL_TAXON_event_data.txt",
                   burnin=0.25, nsamples = 1000)
shifts_n =  credibleShiftSet(edn, 
                             expectedNumberOfShifts=1, 
                             threshold=5, set.limit = 0.95)
# get most likely shift
keynode2 = shifts_n$shiftnodes[[1]]
# get probability 
shifts_n$frequency[1]
r2 = getTipRates(edn)
r2 = data.frame(rate = r2$lambda.avg,
               tip = names(r2$lambda.avg), 
               type = "morphological")

t2 = ggtree(tn, ladderize = F, color = tcols["morphological"], linewidth = 0.3) +
  geom_nodepoint(aes(subset = node %in% keynode2),
                 alpha = 0.5, size = 5)
for (i in c("Ctenotus", "Lerista")) {
  cnode = phytools::findMRCA(tn, tn$tip.label[ grep(i, tn$tip.label) ])
  t2 = t2 + geom_cladelab(node=cnode, label = i, 
                          align=TRUE, 
                          offset = .2, textcolor = tcols["morphological"],
                          barcolor = tcols["morphological"]) 
}
t2 = t2 + xlim(0, 30)

############ 
# clusters
#############

edl = getEventData(tl, 
                   eventdata = "data/diversification/CLUSTERS_event_data.txt",
                   burnin=0.25, nsamples = 1000)
shifts_l =  credibleShiftSet(edl, 
                             expectedNumberOfShifts=1, 
                             threshold=5, set.limit = 0.95)
# get most likely shift
keynode3 = shifts_l$shiftnodes[[1]]
# get probability 
shifts_l$frequency[1]
r3 = getTipRates(edl)
r3 = data.frame(rate = r3$lambda.avg,
                tip = names(r3$lambda.avg), 
                type = "incipient")

t3 = ggtree(tl, ladderize = F, color = tcols["incipient"], linewidth = 0.3) +
  geom_nodepoint(aes(subset = node %in% keynode3),
                 alpha = 0.5, size = 5)
for (i in c("Ctenotus", "Lerista")) {
  cnode = phytools::findMRCA(tl, tl$tip.label[ grep(i, tl$tip.label) ])
  t3 = t3 + geom_cladelab(node=cnode, label = i, 
                          align=TRUE, 
                          offset = .2, textcolor = tcols["incipient"],
                          barcolor = tcols["incipient"]) 
}
t3 = t3 + xlim(0, 30)


############ 
# threshold
#############

edt = getEventData(tt, 
                   eventdata = "data/diversification/THRESHOLD_TAXON_event_data.txt",
                   burnin=0.25, nsamples = 1000)
shifts_t =  credibleShiftSet(edt, 
                             expectedNumberOfShifts=1, 
                             threshold=5, set.limit = 0.95)
# get most likely shift
keynode4 = shifts_t$shiftnodes[[1]]
# get probability 
shifts_t$frequency[1]
r4 = getTipRates(edt)
r4 = data.frame(rate = r4$lambda.avg,
                tip = names(r4$lambda.avg), 
                type = "threshold")

t4 = ggtree(tt, ladderize = F, color = tcols["threshold"], linewidth = 0.3) +
  geom_nodepoint(aes(subset = node %in% keynode4),
                 alpha = 0.5, size = 5)
for (i in c("Ctenotus", "Lerista")) {
  cnode = phytools::findMRCA(tt, tt$tip.label[ grep(i, tt$tip.label) ])
  t4 = t4 + geom_cladelab(node=cnode, label = i, 
                          align=TRUE, 
                          offset = .2, textcolor = tcols["threshold"],
                          barcolor = tcols["threshold"]) 
}
t4 = t4 + xlim(0, 30)

##########
# put it together
##########

rates = rbind(r1, r2, r3, r4)
rates$genus = gsub("_\\S+", "", rates$tip)
rates$genus2 = ifelse(rates$genus %in% c("Ctenotus", "Lerista"), rates$genus, "other")

rates$type2 = factor(rates$type, levels = c("morphological", "operational",  "threshold", "incipient"))


# https://wilkelab.org/ggtext/
d = ggplot(rates, aes(genus2, rate, colour = type2)) +
  geom_boxplot(outlier.colour = NULL) + 
  xlab("") + ylab("speciation rate") + scale_y_log10() +
  theme_classic() +
  scale_colour_manual(values = tcols) +
  theme(legend.title = element_blank(), 
        legend.position = "none")

ab = a + (t2 / t1 / t3 / t4) + d + plot_annotation(tag_levels = 'A') +
   plot_layout(widths = c(8, 5, 8))
cowplot::save_plot("figures/diversification.png", 
                   ab, base_height = 3.5, base_width = 10)
cowplot::save_plot("figures/clads_rates.png", 
                   b, base_height = 2.5, base_width = 8.5)

pdf("figures/operational_css.pdf", height = 4, width = 3.5)
plot.credibleshiftset(shifts_c)
dev.off()

pdf("figures/nominal_css.pdf", height = 4, width = 3.5)
plot.credibleshiftset(shifts_n)
dev.off()

pdf("figures/incipient_css.pdf", height = 4, width = 3.5)
plot.credibleshiftset(shifts_l)
dev.off()


pdf("figures/threshold_css.pdf", height = 4, width = 3.5)
plot.credibleshiftset(shifts_t)
dev.off()

