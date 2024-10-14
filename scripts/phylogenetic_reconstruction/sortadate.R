x1 = read.table("~/Desktop/sortadate/sortadate2.txt")
names(x1) = c("locus", "bp")
x1$locus = gsub(".SH.tre", "", x1$locus)
x2 = read.table("~/Desktop/sortadate/sortadate1.txt")
names(x2) = c("locus", "var", "treelength")
x2$locus = gsub(".SH.tre", "", x2$locus)
x3 = read.csv("~/Desktop/locus_data.csv")

library(tidyverse)
x = full_join(x1, x2)
x = full_join(x, x3)

xx = x[which(x$missingness > 0.95), ]
xx = xx[which(xx$bp > 0.145), ]
xx = xx[which(xx$var < 0.0005), ]
xx = xx[which(xx$treelength < 5), ]

d = data.frame(locus = sample(xx$locus, 25))
write.csv(d, "~/Desktop/sqcl_loci.csv",
          row.names = F)
