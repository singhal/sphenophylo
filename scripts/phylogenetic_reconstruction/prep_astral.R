library(ape)
library(phangorn)

miss = 0.3
AHE = TRUE
tol = 1e-5
collapse = 80

loc = read.csv("/home/babs/spheno/phylogeny/locus_data.csv")
loc = loc[which(loc$missingness > 0.3), ]
if (AHE) {
   loc = loc[grep("AHE", loc$locus), ]
   }
loci = loc$locus

if (AHE) {
   out = paste0("/home/babs/spheno/phylogeny/astral/genetrees_miss", miss, "_tol", tol, "_collapse", collapse, ".AHE.trees")
   } else {
   out = paste0("/home/babs/spheno/phylogeny/astral/genetrees_miss", miss, "_tol", tol, "_collapse", collapse, ".trees")
}	 

trees = vector("list", length(loci))

for (ix in 1:length(loci)) {
    loc = loci[ix]
    t = paste0("/home/babs/spheno/phylogeny/trim_300_0.03_0.3/", loc, ".fasta.aln.treefile")
    if (file.exists(t)) {
       t = read.tree(t)
       t = di2multi(t, tol=tol)
       t$node.label = as.numeric(t$node.label)
       t = pruneTree(t, collapse)

       trees[[ix]] = t
       }
    }
trees = trees[!is.na(trees)]
class(trees) = "multiPhylo"
write.tree(trees, out)