###############################################
# making heatmaps for spot-based CNV analysis #
###############################################

#args = commandArgs(trailingOnly = TRUE)

dir = "/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2/roi_based_cnvcalling" # args[1] #
ext = "LN0438_MAAFHY1_R1" # args[3] # 
cmat = paste0(dir, "/", ext, "/NEW_SPOT_04_CHAR_MAT_2_from_NEW_SPOT_06_cleaned_cnvs_per_bc_0.9.qs") # args[2] # 



# load libraries
library(qs)
library(circlize)
library(GenomicRanges)
library(EnrichedHeatmap)
library(tidyverse)
library(ComplexHeatmap)
library(Matrix)
library(RColorBrewer)


# load spatial domain annotations per barcode, replace "-" by "."
spd = qread("/g/saka/Tatjana/data/all_vis_bs_and_spatial_domain.qs", nthreads = 32)
# only get patient sample needed right now
# when doing roi-based analysis: extract just patient ID from "ext", to subset spatial domains by patient
pat_id = sub("_.*", "", ext)


spd = spd[[pat_id]]
rownames(spd) = gsub("-", ".", rownames(spd))

# load character matrix
char_mat = qread(cmat, 
                 nthreads = 32)

## TEST: keep only chr1 to chr4 in heatmap since the heatmap plots look quite weird
# Keep only rows with chr1, chr2, chr3, or chr4
rows_to_keep <- grepl("chr[1-4]", rownames(char_mat))

# Subset the matrix
char_mat = char_mat[rows_to_keep, ]

# ROI - based analysis:
# subset total barcodes & spatial domains by the bc of specific roi
spd = spd[colnames(char_mat), , drop = FALSE]


domcols = c(
  "follicular_fdc.lo" = "#E6194B",
  "follicular_fdc.hi" = "#3CB44B",
  "interfollicular_malig.lo" = "#FFE119",
  "interfollicular_malig.hi" = "#911EB4",
  "diffuse_fdc.lo" = "#F58231",
  "diffuse_fdc.hi" = "#42D4F4",
  "other" = "#F032E6")


# to color the heatmap by spatial domains, get unique domains per sample
doms_sample = unique(spd[,1])

# subset "domcols" to only keep only colors that are needed for that very sample
# since no sample has all spatial domains present
domcols_sub = domcols[doms_sample]


# transpose so rows = Visium spots
char_mat_t = t(char_mat)
#rm(char_mat)
#gc()

# extract chromosome from colnames of char mat, get unique chroms
chroms = sub(":.*", "", colnames(char_mat_t))
chroms_u = unique(chroms)
chroms_o = order(chroms)

# IMPORTANT: set factor lvls manually, otherwise chromosomes are sorted alphanumerically when plotting (chr1, chr10, chr11, ...)
# need to do to display "char_mat_t" correctly in heatmap
chrom_lvls = c(paste0("chr", 1:22), "chrX", "chrY")
chroms_fact = factor(chroms, levels = chrom_lvls)

# get bin start pos from colnames
starts = as.numeric(sub(".*:([0-9]+)-([0-9]+)", "\\1", colnames(char_mat_t)))

# make granges for bins (for creating ideogram)
gr = GRanges(seqnames = chroms,
             ranges = IRanges(start = starts, 
                              end = starts + 9999))  # assuming 10kb bins

# set cnv colors
cnv_colors = c("amp" = "red", "del" = "blue", "NA" = "white") # very sloppy, NA vals should be NA and not string NA - quick workaround, fix later
heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12),   # legend title
    labels_gp = gpar(fontsize = 10)   # legend labels
  )

# make row annotation for spatial domains
row_anno = rowAnnotation(
  SpatialDomain = spd[,1],
  col = list(SpatialDomain = domcols_sub),
  annotation_name_side = "top",
  show_annotation_name = TRUE)

# make annot for idogram: each chrom gets bar which spans its bins
ideo_anno = HeatmapAnnotation(
  ideogram = anno_barplot(
    rep(1, ncol(char_mat_t)),   
    gp = gpar(fill = "#BEBEBE", col = NA),
    border = TRUE,
    axis = FALSE,
    annotation_name = NULL,
    column_title_gp = gpar(fontsize = 8)
  ),
  annotation_name = NULL,
  height = unit(0.5, "cm")
)

#png(paste0(dir, "/", ext, "/NEW_SPOT_05_heatmap_from_NEW_SPOT_03_75thqn_cnvs_20251117.png"), width = 8700, height = 2300, res = 400)

# now, make heatmap with column split to get individual "boxes" for chromosomes
htmp = Heatmap(
  char_mat_t,
  name = "CNV",
  col = cnv_colors,
  na_col = "white", ##
  border = TRUE,
  row_title = "barcoded spots",
  row_title_gp = gpar(fontsize = 17),
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean", # added 17.11.25
  clustering_method_rows = "complete", # added 17.11.25
  show_row_dend = TRUE, # added 17.11.25
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_split = chroms_fact,
  column_title = chroms_u,
  column_title_gp = gpar(fontsize = 15),
  bottom_annotation = ideo_anno,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 16), 
    labels_gp = gpar(fontsize = 12)
    ), 
  left_annotation = row_anno,
  column_order = colnames(char_mat_t)
  # width = rep(unit(1, "cm"), length(chroms_u)) # heatmap not being plotted when I do this, look for it later
)

draw(
  htmp, 
  heatmap_legend_side = "bottom", 
  column_title = paste(ext, "- filtered CNV profile (75th quantile threshold)"), 
  column_title_gp = gpar(fontsize = 23)
  )

#dev.off()
