library(phangorn)
library(dplyr)
library(ape)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Ctenotus_speciation/")
d = read.csv("metadata/sphenomorphine_all_v4.csv")

d1 = d[which(d$SQCL == TRUE), ]
d1 = d1[which(d1$CYTB == TRUE), ]
d1 = d1[which(d1$ddRAD == TRUE),]
d2 = d1 %>% select(SAMPLE_ID, PROBES, CLUSTER_SUBGENUS)
d3 = d2 %>% group_by(CLUSTER_SUBGENUS) %>% 
  arrange(PROBES) %>% slice_tail(n  = 1)

grs = unique(d[complete.cases(d$CLUSTER_SUBGENUS), "CLUSTER_SUBGENUS"])

tx = read.tree("phylogeny/SqCL/concat/concat_ind0.05_loci0.3_all_n5277.phy.treefile.tre")
outs = d[which(d$INGROUP == FALSE), "SAMPLE_ID"]
outs = c(outs, "spheno_unk1", "spheno_unk3")
tx = root(tx, outs[1], resolve.root = T)
tx = drop.tip(tx, outs)
drop = c("UMMZ_244315_ct_quat", "CUMV_14452_Le_bipe")
tx = drop.tip(tx, drop)
tx1 = keep.tip(tx, d3$SAMPLE_ID)
tx1$tip.label = pull(d3[ match(tx1$tip.label, d3$SAMPLE_ID), "CLUSTER_SUBGENUS"])

for (i in 1:length(grs)) {
  gr = grs[i]
  grix = which(tx1$tip.label == gr)
  
  nn = phangorn::Ancestors(tx1, grix)
  
  tips = lapply(nn, function(x) {tx1$tip.label[ phangorn::Descendants(tx1, x, type = "tips")[[1]] ]})
  tips = unique(unlist(tips))
  tips = tips[tips != gr][1:3]
  tips = pull(d3[d3$CLUSTER_SUBGENUS %in% tips, "SAMPLE_ID"])
  tips2 = paste(tips, collapse = ",")
  
  outtips = paste0(tips, ".outgroup")
  tipx = c(d[which(d$CLUSTER_SUBGENUS == gr), "SAMPLE_ID"], tips)
  contree = keep.tip(tx, tipx[ tipx %in% tx$tip.label ])
  contree$tip.label[match(tips, contree$tip.label)] = outtips
  
  if (Ntip(contree) > 3) {
    conout = paste0("~/Desktop/", gr, ".tre")
    write.tree(contree, conout)
    
  } else {
    cat(gr, "\n")
  }
  call = paste0("python /media/babs/brains3/spheno/scripts/create_alignment_reference_bias_inclusive_RAM.py ",
                  "--miss 0.6 --level CLUSTER_SUBGENUS --name ", grs[i],
                  " --outsp ", tips2, "\n")
  call1 = paste0("python /media/babs/brains3/spheno/scripts/create_alignment_reference_bias_inclusive_RAM.py ",
                   "--miss 0.6 --level CLUSTER_SUBGENUS --name ", grs[i],
                   " --outsp ", tips2, " --mito\n")
  
  cat(call)
  cat(call1)
  
}


