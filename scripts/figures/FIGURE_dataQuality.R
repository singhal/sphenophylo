rm(list = ls())
library(ggplot2)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
d = read.csv("../Sphenomorphine_Phylogeny/data/tip_data.csv")
t = read.tree("../Sphenomorphine_speciation/phylogeny/SqCL/concat/concat_ind0.05_loci0.3_all_n5277.phy.treefile.tre")
x = read.csv("../Sphenomorphine_speciation/metadata/sphenomorphine_all_v10.csv")

d1 = d[d$ind %in% t$tip.label, ]
# all the missing ones are outgroups
miss = t$tip.label[ ! t$tip.label %in% d$ind ]

d1$prev_published = FALSE

# previously published data
x = read.csv("metadata/previously_published_SqCL/Title_etal.csv")
d1[(d1$ind %in% x$SAMPLE), "prev_published"] = TRUE

y = read.table("metadata/previously_published_SqCL/Grundler_etal.txt",
               sep = "\t", header = T)
for (i in 1:nrow(y)) {
  m = grep(y[i, "sample"], d1$ind)
  if (length(m) == 1) {
    d1[m, "prev_published"] = TRUE
  } else {
    cat(y[i, "sample"], "\n")
  }
}
# outgroups that somehwo went uncaught - 
d1[d1$ind == "FMNH258830", "prev_published"] = TRUE
d1[d1$ind == "FMNH258750", "prev_published"] = TRUE
d1[d1$ind == "spheno_unk3", "prev_published"] = TRUE
table(d1$prev_published)
d2 = d1 %>% filter(prev_published == FALSE)
d2$num_sampled = ifelse(d2$num_UCE_contigs < 100, 394, 5462)
d2$per_sampled = d2$num_total / d2$num_sampled

a = ggplot(d2, aes(per_sampled)) + 
  geom_histogram() + xlab("% of loci sampled") +
  theme_classic()
b = ggplot(d2, aes(as.numeric(avg_length))) + 
  geom_histogram() + xlab("avg. locus length (bp)") +
  theme_classic()
c = ggplot(d2, aes(as.numeric(avg_cov))) + 
  geom_histogram() + xlab("avg. coverage") +
  theme_classic()
d = ggplot(d2, aes(as.numeric(heterozygosity))) + 
  geom_histogram() + xlab("heterozygosity") +
  theme_classic()

qual = cowplot::plot_grid(plotlist = list(a, b, c, d), nrow = 2, 
                   labels = c("A", "B", "C", "D"))
cowplot::save_plot("figures/data_quality.png", qual, 
                   base_height = 4, base_width = 6)
