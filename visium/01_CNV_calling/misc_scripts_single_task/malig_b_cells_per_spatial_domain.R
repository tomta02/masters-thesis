############################################
# plot malignant B cells per spatial domain
############################################

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


# add plot title name!!!.
plotlist = list()
for (i in (seq_along(visiumlist))) {
  v = visiumlist[[i]]
  n = names(visiumlist[i])
  
  cell_abun = v@assays$cell2loc_q05$counts
  cell_abun = as.data.frame(t(cell_abun))
  
  spd = v@meta.data$prelim_growth_2
  
  cell_abun = cell_abun %>%
    mutate(
      row_sums = rowSums(across(everything())),
      malig_sum = rowSums(across(contains("malig"))),
      malig_frac = malig_sum/row_sums
    )
  
  plot_df = data.frame(
    malig_frac = cell_abun$malig_frac, 
    spatial_domain = spd
  )
  
  
  
  p = ggplot(plot_df, aes(x = spatial_domain, y = malig_frac, fill = spatial_domain)) +
    geom_violin(trim = FALSE, color = "black") +
    stat_summary(
      fun = mean, 
      geom = "point", 
      shape = 15, 
      size = 2, 
      color = "black",
      show.legend = FALSE
    ) +
    stat_summary(
      fun.data = mean_sdl, 
      fun.args = list(mult = 1),  # 1 standard deviation
      geom = "errorbar", 
      width = 0.2, 
      color = "black"
    ) +
    #geom_jitter(size = 1, alpha = 0.3) +
    # geom_jitter(shape = 21, size = 1, alpha = 0.3, color = NA) + 
    scale_fill_manual(values = domcols) +
    theme_classic() +
    ggtitle(n) +
    xlab("Spatial Domain") +
    ylab("malignant B cell fraction") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  plotlist[[i]] = p
}


for (i in (seq_along(plotlist))) {
  n = names(visiumlist[i])
  ggsave(
    paste0("/g/saka/Tatjana/data/05_plots/malig_b_cell_fraction_per_roi/pdf/", n, "_malignant_B_cell_frac.pdf"),
    plotlist[[i]], 
    width = 7,
    height = 4,
    units = "in",
    dpi = 300
  )
}





