###############################################
# making heatmaps for spot-based CNV analysis #
###############################################

args = commandArgs(trailingOnly = TRUE)

dir = args[1]
cmat = args[2]
ext = args[3]


# load libraries
library(qs)
library(circlize)
library(GenomicRanges)
library(EnrichedHeatmap)
library(tidyverse)
library(ComplexHeatmap)
library(Matrix)


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

# ROI - based analysis:
# subset all barcodes & spatial domains of that patient by the specific ROI
spd = spd[colnames(char_mat), , drop = FALSE]

# to later color the heatmap by spatial domains, get all unique ones
doms = unique(spd[,1])

# generate random colors for now, coordinate color coding w Kristy later
domcols = setNames(rand_color(length(doms)), doms)
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

# make row annotation for spatial domains
row_anno = rowAnnotation(
  SpatialDomain = spd[,1],
  col = list(SpatialDomain = domcols),
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

png(paste0(dir, "/", ext, "/", "SPOT_05_heatmap_clustrows_TRUE.png"), width = 8500, height = 2000, res = 300)

# now, make heatmap with column split to get individual "boxes" for chromosomes
htmp = Heatmap(
  char_mat_t,
  name = "CNV",
  col = cnv_colors,
  na_col = "white", ##
  border = TRUE,
  row_title = "Visium spots",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_split = chroms_fact,
  column_title = chroms_u,
  bottom_annotation = ideo_anno,
  left_annotation = row_anno,
  #column_order = colnames(char_mat_t)
  # width = rep(unit(1, "cm"), length(chroms_u)) # heatmap not being plotted when I do this, look for it later
)

draw(
  htmp, 
  heatmap_legend_side = "bottom", 
  column_title = ext, 
  column_title_gp = gpar(fontsize = 15)
  )

dev.off()
