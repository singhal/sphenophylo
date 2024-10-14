do_ancestors = function(phy) {
  paths = ape::nodepath(phy)
  function(tip) {
    return (rev(head(paths[[tip]], -1)))
  }
}

do_tips = function(phy) {
  paths = ape::nodepath(phy)
  function(ancestor) {
    which(sapply(paths, function(p) ancestor %in% p))
  }
}

do_brlens = function(phy) {
  len = numeric(ape::Ntip(phy) + phy$Nnode)
  len[phy$edge[,2]] = phy$edge.length
  function(node) {
    len[node]
  }
}

DR = function(phy) {
  brlens = do_brlens(phy)
  tips = do_tips(phy)
  ancestors = do_ancestors(phy)
  
  x = numeric(ape::Ntip(phy))
  for (i in 1:ape::Ntip(phy)) {
    p = c(i, ancestors(i))
    for (j in seq_along(p))
      x[i] = x[i] + brlens(p[j]) / 2^(j-1)
  }
  drrates = 1 / x
  names(drrates) = phy$tip.label
  return(drrates)
}
