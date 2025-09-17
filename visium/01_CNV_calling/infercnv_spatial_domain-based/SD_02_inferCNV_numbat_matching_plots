###############################################
# CNV matching numbat and infercnv - plots    
###############################################

args = commandArgs(trailingOnly = TRUE)

nb = args[1]
infercnv = args[2]
matching = args[3]
outpath = args[4]

# load libraries
library(qs)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(GenomicRanges)
library(VennDiagram)
library(gtrellis)
library(circlize)

# load data:
# cnvs called from numbat in the single cell data
# cnvs called from infercnv in visium data
# cnvs that were shared/matching between the 2
numbat = read.csv(nb)
infcnv = read.csv(infercnv)
shared = read.csv(matching)

# modify data... 
# gtrellis expects the integer nr. of the chromosome to be preceded by "chr"
# fix in other script later
shared$chr <- paste0("chr", shared$chr)
numbat$CHROM <- paste0("chr", numbat$CHROM)
sicnv$chr <- paste0("chr", sicnv$chr)

# small helper function to lmake GRanges objects
make_granges = function(df, chr_col, start_col, end_col, state_col) {
  GRanges(
    seqnames = df[[chr_col]],
    ranges = IRanges(start = df[[start_col]], end = df[[end_col]]),
    cnv_state = df[[state_col]]
  )
}

gr_infcnv = make_granges(infcnv, "chr", "start", "end", "cnv_state_clean")
gr_numbat = make_granges(numbat, "CHROM", "seg_start", "seg_end", "cnv_state_post")
gr_shared = make_granges(shared, "chr", "start_match", "end_match", "state")

# make gtrellis plot
pdf(outpath, width = 15, height = 4)

# for some reason, this does not work
col_fun = c("amp" = "red", "del" = "blue")


gtrellis_layout(
  add_name_track = TRUE,
  n_track = 3, 
  equal_width = TRUE,
  title = "CNVs inferred using inferCNV, numbat, and overlapping CNVs", 
  track_ylab = c("inferCNV", "numbat", "overlap"), 
  track_axis = FALSE,
  xlab = "")

add_rect_track(
  gr_sicnv,
  track = 2,
  h1 = 1, 
  h2 = 0,
  gp = gpar(col = NA, fill = col_fun[mcols(gr_sicnv)$cnv_state])
)

add_rect_track(
  gr_numbat,
  track = 3,
  h1 = 1, 
  h2 = 0,
  gp = gpar(col = NA, fill = col_fun[mcols(gr_numbat)$cnv_state])
)

add_rect_track(
  gr_shared,
  track = 4,
  h1 = 1, 
  h2 = 0,
  gp = gpar(col = NA, fill = col_fun[mcols(gr_shared)$cnv_state])
)

dev.off()
