# plotting CNV information on top of spatial transcriptomics data
##############################################################

#args = commandArgs(trailingOnly = TRUE)

dir = "/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2/roi_based_cnvcalling" # args[1] #   
#cnv = args[2] # paste0(dir, "/LN0027_1VL5B2_R4/NEW_SPOT_03_cnvs_cleaned_75th_quantile_binarized.csv") #
ext = "LN0027_1VL5B2_R4" # args[4] # 
cnv = (paste0(dir, "/", ext, "/NEW_SPOT_03_cnvs_cleaned_75th_quantile_binarized.csv"))      # args[3]
#outdir = "/g/saka/Tatjana/data/01_CNV_analysis/03_plots/cnv_content_in_space/"

# libraries
library(Seurat)
library(tidyverse)
library(qs)
library(GenomicRanges)

# load visium objects
visiumlist = qread("/g/saka/Kristy/projects/composite/analysis/rdata/02_13_visium_flt.qs", nthreads = 32)

# load cnv info
cnv = read.csv(cnv)
# important for later matching - cnv barcode names (i.e. colnames) contain "." a
# nd the visium barcode names contain "-". replace
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



# # now plot spatialplot
# p = SpatialFeaturePlot(
#   v, 
#   "CNV_bp_full",
#   pt.size.factor = 2.5,
#   image.scale = "hires",
#   image.alpha = 0.2
#   ) +
#   ggtitle(paste0("CNV-derived genomic content: ", ext)) +
#   labs(fill = "CNV bp (scaled)") +
#   theme(aspect.ratio = 1)
# 
# ggsave(
#   paste0(outdir, ext, "_cnv_bp_in_space.pdf"),
#   p,
#   dpi = 300)


# # maybe also try out with binarized CNVs?
# cnv_bin = read.csv(cnv_bina)
# colnames(cnv_bin) = gsub("\\.", "-", colnames(cnv_bin))
# 
# # make simple summary of count chunks
# cnv_res_bin = data.frame(
#   bc = names(cnv_bin)[-c(1:6)],
#   cnv_frac = colSums(!is.na(cnv_bin[ , -c(1:6)]))
# )
# 
# # get max value
# m2 = max(cnv_res_bin$cnv_frac)
# 
# cnv_res_bin$scaled = (cnv_res_bin$cnv_frac/m2)
# #cnv_res_bin$cnv_frac = NULL
# 
# # add CNV information as metadata to visium object
# v[["CNV_bin"]] = cnv_res_bin[colnames(v), "scaled"]
# 
# # now plot spatialplot
# p2 = SpatialFeaturePlot(
#   v,
#   "CNV_bin",
#   pt.size.factor = 2.5,
#   image.scale = "hires",
#   image.alpha = 0.2
# ) +
#   ggtitle(paste0("bin. CNV-derived genomic content: ", ext)) +
#   labs(fill = "CNV bp (scaled)") +
#   theme(aspect.ratio = 1)
# 
# ggsave(
#   paste0(outdir, ext, "_perspotbin_cnv_bp_in_space.pdf"),
#   p2,
#   dpi = 300)
# 
