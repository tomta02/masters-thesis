##########################################
# plotting cell2loc cell type abundances
##########################################

# load libraries
library(Seurat)
library(tidyverse)
library(qs)
library(patchwork)

# dir for saving figs later
fig = "/g/saka/Tatjana/data/05_plots/c2l_celltype_abundance/NEW/"

# load visium objects
visiumlist = qread("/g/saka/Kristy/projects/composite/analysis/rdata/02_13_visium_flt.qs", nthreads = 32)

#test = visiumlist[[1]]
# looking at cell2loc celltypes to figure out which 
# groups of cells to group together for the 
# celltype comparison to codex
#c = test@assays$cell2loc_q05$counts
#c = t(c)

# make featurs to plot in spatialfeatureplot,
# for looping over
features = c("B_frac", "T_frac", "DC_frac", "FRC_frac")
 
pltlist = lapply(visiumlist, function(v) {
  ca = v@assays$cell2loc_q05$counts
  ca = as.data.frame(t(ca))
  
  ca = ca %>%
    mutate(
      all_sum = rowSums(across(everything())),
      B_sum = rowSums(across(contains("B-"))) +
        rowSums(across(contains("malig"))),
      T_sum = rowSums(across(contains("T-"))),
      DC_sum = rowSums(across(contains("DC"))),
      FRC_sum = rowSums(across(contains("VSMC"))),
      
      # now caluclate ratios
      B_frac = (B_sum/all_sum),
      T_frac = (T_sum/all_sum),
      DC_frac = (DC_sum/all_sum),
      FRC_frac = (FRC_sum/all_sum)
      
    )
  
  #add to seurat metadata for poltting
  v = AddMetaData(v, ca[, features])
  
  # makee spatialfeatureplots
  p = lapply(features, function(f) {
    SpatialFeaturePlot(
      v,
      features = f,
      pt.size.factor = 2.5,
      image.scale = "hires",
      image.alpha = 0.2) + 
      theme(aspect.ratio = 1)
  })
  p = wrap_plots(p, nrow = 1)
})

# save plots
for (i in seq_along(pltlist)) {
  p = pltlist[[i]]
  ggsave(
    paste0(fig, i, "_c2l_celltype_dist.png"),
    p,
    height = 6,
    width = 11,
    dpi = 300
  )
}

