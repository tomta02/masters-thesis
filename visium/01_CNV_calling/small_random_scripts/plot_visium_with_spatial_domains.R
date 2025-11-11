#########################################
# PLOT VISIUM ROIS AND SPATIAL DOMAINS
#########################################

# load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(qs)
library(RColorBrewer)

# load visium objects
visiumlist = qread("/g/saka/Kristy/projects/composite/analysis/rdata/02_13_visium_flt.qs", nthreads = 32)


# set colors for spatial domains - set manually instead, did not likee color palettes from colorbrewer
domcols = c(
  "follicular_fdc.lo" = "#E6194B",
  "follicular_fdc.hi" = "#3CB44B",
  "interfollicular_malig.lo" = "#FFE119",
  "interfollicular_malig.hi" = "#911EB4",
  "diffuse_fdc.lo" = "#F58231",
  "diffuse_fdc.hi" = "#42D4F4",
  "other" = "#F032E6")

# make separate "legend" plot for spatial domains and associated colors
domcols_df = data.frame(
  spd = names(domcols),
  x = 1,
  y = seq_along(domcols)
)

legend_p = ggplot(domcols_df, aes(x = x, y = y, fill = spd)) +
  geom_tile(width = 0.8, height = 0.8) +
  scale_fill_manual(values = domcols) + 
  geom_text(aes(label = spd, y = y), x = 1.8, hjust = 0, size = 5) +
  theme_void() +
  theme(legend.position = "none") +
  xlim(0.5, 3)  

ggsave(
  "/g/saka/Tatjana/data/05_plots/spatial_doms_legend.png",
  legend_p,
  dpi = 400)



# to keep list of plots instead of saving directyl
plotlist = list()

for (i in seq_along(visiumlist)) {
  roi = visiumlist[[i]]
  roiname = names(visiumlist[i])
  
  p = SpatialPlot(
    roi, 
    group.by = "prelim_growth_2",
    cols = domcols,
    image.scale = "hires",
    image.alpha = 0.1,
    label.size = 10
    ) + ggtitle(roiname) + 
     theme(
      plot.title = element_text(hjust = 0.5, size = 22),
      #legend.text = element_text(size = 12),
      #legend.title = element_text(size = 14),
      #legend.key.size = unit(1.2, "cm")
      legend.position = "none",
      aspect.ratio = 1)
  
  plotlist[[i]] = p
}


for (i in seq_along(plotlist)) {
  ggsave(
    paste0("/g/saka/Tatjana/data/05_plots/", i, "_spatial_domain_plot_nolegend.png"),
    plotlist[[i]],
    width = 5,
    height = 5,
    dpi = 400)
}


