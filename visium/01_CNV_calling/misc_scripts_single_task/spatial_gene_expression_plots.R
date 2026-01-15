####################################
# Spatial plots of gene expression
####################################

# load libraries
library(Seurat)
library(tidyverse)
library(qs)

# dir for saving figs
fig = "/g/saka/Tatjana/data/05_plots/visium_gene_expression/all_4/"

# load visium objects
visiumlist = qread("/g/saka/Kristy/projects/composite/analysis/rdata/02_13_visium_flt.qs", nthreads = 32)

# current default assay is cell2loc output, switch to SCT
visiumlist = lapply(visiumlist, function(v) {
  DefaultAssay(v) = "SCT"
  v
})

#-------------------------------------------------------------------------------
# check the expression of following marker genes for comparison w codex:
# B-cells : PAX5, CD20       | gene name: PAX5, MS4A1
# T-cells: CD3               | gene name: CD3E, CD3D, CD3G, [CD3Z not in tissue]
# FRCs: CD90                 | gene name: THY1
# Dendritic cells: CD11c     | gene name: ITGAX


test = visiumlist[[4]]

# create module scores for combined gene epression 
# PAX5 and CD20 as proxy for B cells, all genes making up T cell receptor
g = list(
  Bcells = c("PAX5", "MS4A1"),
  Tcells = c("CD3E", "CD3D", "CD3G"),
  FRC = c("THY1"),
  DC = c("ITGAX")
)

# add modulescores to visium objects
visiumlist = lapply(visiumlist, function(v) {
  AddModuleScore(
    v,
    features = g,
    name = names(g)
  )
})


# make patchwork fig: expression of each of the 4 module scores, for 
# 1 sample at a time
# the module scores are callled the sam ein each sample: 
# "Bcells1", "Tcells2", "FRC3", "DC4"
scores = c("Bcells1", "Tcells2", "FRC3", "DC4")

pl = list()

for (i in (seq_along(visiumlist))) {
  v = visiumlist[[i]]
  
  p = lapply(scores, function(x) {
    SpatialFeaturePlot(
      v, 
      features = x,
      pt.size.factor = 2.5,
      image.scale = "hires",
      image.alpha = 0.2) +
      theme(aspect.ratio = 1)
  })

  #c = ((p[[1]] | p[[2]]) / (p[[3]] | p[[4]]) +
    #plot_layout(widths = c(1, 1), heights = c(1, 1)))
  combined = patchwork::wrap_plots(p, ncol = 2, nrow = 2)
  
  pl[[i]] = combined
}

# save plots
for (i in seq_along(pl)) {
  p = pl[[i]]
  ggsave(
    paste0(fig, i, "_vis_genexpression_ct_abundance_noFRC.png"),
    p,
    height = 8,
    width = 16,
    dpi = 300
  )
}


