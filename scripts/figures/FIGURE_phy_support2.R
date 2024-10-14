rm(list = ls())

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
d = read.csv("../Sphenomorphine_speciation/metadata/sphenomorphine_all_v11.csv")

t1 = read.tree("data/phylogenetic_support/synthetic_tree.all.support.dated.tre")
outs = d[which(d$INGROUP == FALSE), "SAMPLE_ID"]
t2 = drop.tip(t1, outs)

table(complete.cases(t2$node.label))
t2$node.label = as.numeric(gsub("/\\S+", "", t2$node.label))
t2$node.label[!complete.cases(t2$node.label)] = median(t2$node.label, na.rm =T)

cols = RColorBrewer::brewer.pal(11, "Spectral")
pal = colorRampPalette(cols)(101)
# 
# edgecols = rep("black", length(t2$edge.length))
# edgewids = rep(0.1, length(t2$edge.length))
# for (i in 1:length(t2$node.label)) {
#   nodenum = Ntip(t2) + i
#   support = t2$node.label[i]
#   
#   edgenums = which(t2$edge[,1] == nodenum) 
#   # edgecols[edgenums] = pal[round(support) + 1]
#   if (support < 90) {
#     edgecols[edgenums] = "red"
#     edgewids[edgenums] = 1
#   }
# }

pdf("figures/support.pdf", width = 10, height = 10)
par(mar = c(0,0,0,0))
plot.phylo(t2, show.tip.label = FALSE, type = "fan", edge.width = 0.5, edge.color = alpha("gray50", 0.7))
for (i in 1:length(t2$node.label)) {
  nodelabels("", Ntip(t2) + i, pch = 16, cex = 0.7, 
             frame = "none", col = pal[round(t2$node.label[i]) + 1])
}

legend_image <- as.raster(matrix(rev(pal), ncol=1))
op <- par(  ## set and store par
  fig=c(grconvertX(c(-24, -12), from="user", to="ndc"),    ## set figure region
        grconvertY(c(-26, -14), from="user", to="ndc")), 
  mar=c(1, 1, 1, 9.5),                                  ## set margins
  new=TRUE)                                ## set new for overplot w/ next plot

plot.window(c(0, 2), c(0, 100))                               ## ini plot2
rasterImage(legend_image, 0, 0, 1, 100)                       ## the gradient
lbsq <- seq.int(0, 100, l=5)                                  ## seq. for labels
axis(4, at=lbsq, pos=1, labels=F, col=0, col.ticks=1, tck=-.1)  ## axis ticks
mtext(lbsq, 4, 0.01, at=lbsq, las=2, cex=.6)                    ## tick labels
mtext('support', 3, -.125, cex=.6, adj=.1, font=2)          ## legend title

par(op)  ## reset par
dev.off()