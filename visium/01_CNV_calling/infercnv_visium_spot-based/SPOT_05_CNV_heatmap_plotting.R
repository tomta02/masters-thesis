###############################################
# making heatmaps for spot-based CNV analysis #
###############################################

## TODO: 
# add spatial domain annotations to plot
# add ROI name ("ext") as header


#args = commandArgs(trailingOnly = TRUE)

dir = "/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2" #args[1]
cmat = "/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2/LN0438_MAAFHY1_R1/cnv_char_mat_for_heatmap.qs" #args[2]
ext = "LN0438_MAAFHY1_R1" #args[3]

# load libraries
library(qs)
library(circlize)
library(GenomicRanges)
library(EnrichedHeatmap)
library(tidyverse)
library(ComplexHeatmap)


char_mat = qread(cmat, 
                 nthreads = 32)


# transpose so rows = Visium spots
char_mat_t = t(char_mat)

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

png(paste0(dir, "/", ext, "/SPOT_04_TESTHEATMAP.png"), width = 8500, height = 2000, res = 300)


# now, make heatmap with column split to get individual "boxes" for chromosomes
htmp = Heatmap(
  char_mat_t,
  name = "CNV",
  col = cnv_colors,
  na_col = "white", ##
  border = TRUE,
  row_title = "Visium spots",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_split = chroms_fact,
  column_title = chroms_u,
  bottom_annotation = ideo_anno,
  column_order = colnames(char_mat_t)
  # width = rep(unit(1, "cm"), length(chroms_u)) # heatmap not being plotted when I do this, look for it later
)

draw(
  htmp, 
  heatmap_legend_side = "bottom", 
  column_title = ext, 
  column_title_gp = gpar(fontsize = 15)
  )

dev.off()
