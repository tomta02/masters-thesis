############################################
# making numbat track plot
# similar to reference cnv plot in SPOT_06
############################################

args = commandArgs(trailingOnly = TRUE)

numbat = args[1] 
genotypes = args[2] 
pat_id = args[3] 


# load libraries
library(qs)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(GenomicRanges)
library(VennDiagram)
library(gtrellis)
library(circlize)


# load list with numbat data, load genotypes
nb = qread(numbat, nthreads = 32)
gt = qread(genotypes, nthreads = 32)


# extract cnv calls from data
nb_cnv_list = list()

for (name in names(nb)) {
  pat = nb[[name]][["segs_consensus"]]
  nb_cnv_list[[name]] = pat
}

# subset nb_cnv_list by the genotypes that were significant
nb_cnv_list_subs = mapply(function(patient, genot) {
  patient[patient$seg_cons %in% genot, ]
}, nb_cnv_list, gt, SIMPLIFY = FALSE)


# subset numbat cnvs
nb_cnv = as.data.frame(nb_cnv_list_subs[[pat_id]])

# make GenomicRanges object for gtrellis plot
nb_cnv_gr = GRanges(seqnames = nb_cnv$CHROM,
                     ranges = IRanges(
                       start = nb_cnv$seg_start,
                       end = nb_cnv$seg_end),
                     state = nb_cnv$cnv_state_post
)

# "seqnames" from numbat are only integer numbers, but gtrellis expects (chr1, chr2, ... instead of: 1, 2)
# change this for plotting
seqlevelsStyle(nb_cnv_gr) = "UCSC"


# #making gtrellis plot
out_png = paste0("/g/saka/Tatjana/data/01_CNV_analysis/03_plots/", pat_id, "_numbat_cnvs.png")
png(out_png, width = 20, height = 2, units = "in", res = 400)

col_fun = c("amp" = "red", "del" = "blue")

gtrellis_layout(
  species = "hg38",
  add_name_track = TRUE,
  add_ideogram_track = TRUE,
  n_track = 1,
  equal_width = FALSE,
  title = paste0("numbat CNVs  - ", pat_id),
  #track_ylab = c(""),
  track_axis = FALSE,
  xlab = "")

add_rect_track(
  nb_cnv_gr,
  track = 2,
  h1 = 1,
  h2 = 0,
  gp = gpar(col = NA, fill = col_fun[mcols(nb_cnv_gr)$state])
)

dev.off()