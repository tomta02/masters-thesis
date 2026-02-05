#####################################################
# plotting scatterplot malig. B cell frac. vs. CNVs 
#####################################################

args = commandArgs(trailingOnly = TRUE)

dir = args[1] # "/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2/roi_based_cnvcalling" # 
cnvs = args[2] # (paste0(dir, "/", ext, "/NEW_SPOT_03_cnvs_cleaned_75th_quantile_binarized.csv"))  
ext = args[4] # "LN0027_1VL5B2_R4" # 
    

# load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(qs)
library(RColorBrewer)

# load visium objects
visiumlist = qread("/g/saka/Kristy/projects/composite/analysis/rdata/02_13_visium_flt.qs", nthreads = 32)

# load cnv info
cnv = read.csv(cnvs)
# important for later matching - cnv barcode names (i.e. colnames) contain "." 
# and the visium barcode names contain "-". replace
colnames(cnv) = gsub("\\.", "-", colnames(cnv))

# sum df$width for each non-NA entry in each column
cnv_res = data.frame(
  bc = names(cnv)[-c(1:6)],
  width_sum = sapply(cnv[ , -c(1:6)], function(x) sum(cnv$width[!is.na(x)]))
)


# get max value
m = max(cnv_res$width_sum)

cnv_res$scaled = (cnv_res$width_sum/m)
cnv_res$width_sum = NULL


v = visiumlist[[ext]]

# add CNV information as metadata to visium object
v[["CNV_bp_full"]] = cnv_res[colnames(v), "scaled"]


cell_abun = v@assays$cell2loc_q05$counts
cell_abun = as.data.frame(t(cell_abun))

spd = v@meta.data$prelim_growth_2

cell_abun = cell_abun %>%
  mutate(
    row_sums = rowSums(across(everything())),
    malig_sum = rowSums(across(contains("malig"))),
    malig_frac = malig_sum/row_sums
  )

v[["maligBfrac"]] = cell_abun[colnames(v), "malig_frac"]


# set colors for spatial domains 
domcols = c(
  "follicular_fdc.lo" = "#E6194B",
  "follicular_fdc.hi" = "#3CB44B",
  "interfollicular_malig.lo" = "#FFE119",
  "interfollicular_malig.hi" = "#911EB4",
  "diffuse_fdc.lo" = "#F58231",
  "diffuse_fdc.hi" = "#42D4F4",
  "other" = "#F032E6")

# calculate regression between two variables
fit = lm(CNV_bp_full ~ maligBfrac, data = v@meta.data)
r2  = summary(fit)$r.squared

# make scatterplot
p = ggplot(v@meta.data, aes(x = maligBfrac, y = CNV_bp_full, color = prelim_growth_2)) +
  geom_point(alpha = 0.4, size = 2)+
  scale_color_manual(values = domcols) +
  theme_minimal() +                                 
  labs(
    x = "malignant B cell fraction",
    y = "relative CNV burden",
    title = ext,
  ) +
  geom_smooth(method = "lm", color = "black", se = FALSE, linewidth = 0.25) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("R² = ", round(r2, 3)),
    hjust = 1.1, vjust = 1.5,
    size = 4,
    color = "black"
  ) +
  theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
)


ggsave(
  paste0("/g/saka/Tatjana/data/05_plots/scatterplots_malig_b_vs_cnv_burden_per_spot_with_R2/", ext, "_scatter.png"),
  p,
  width = 8,
  height = 4,
  dpi = 300
)


