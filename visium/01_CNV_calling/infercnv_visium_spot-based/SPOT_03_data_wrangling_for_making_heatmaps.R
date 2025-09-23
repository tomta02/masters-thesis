#########################################
# data wrangling for heatmap generation #
#########################################

args = commandArgs(trailingOnly = TRUE)

dir = args[1]
cnvs = args[2]
ext = args[3]

########################################################################################
# generate character matrix of cnv data to plot cnvs that were uncovered per visium spot
# as opposed to per "spatial domain".
# using ComplexHeatmap library described here:
# https://jokergoo.github.io/2020/10/29/make-genome-scale-heatmap/

# load libraries
library(qs)
library(gtrellis)
library(circlize)
library(GenomicRanges)
library(EnrichedHeatmap)
library(tidyverse)
library(ComplexHeatmap)

cnv = read.csv(cnvs)

# the "normalize_genomic_signals_to_bins" function does not work with a full character matrix (only numeric matrices)
# -> need to supply our data as individual vectors. To do this more easily, split our cnv data per visium spot barcode
# and save thesee "slices" in list
# there's prob a way better way to do this, think about
base_cols = cnv[, 1:6] # are metadata cols, exclude

# get all the patient columns
pat_id = sub("_.*", "", ext)
bc_cols = grep(pat_id, names(cnv), value = TRUE)

bc_sublists = map(bc_cols, function(col_name) {
  cbind(base_cols, cnv[[col_name]]) %>%
    setNames(c(names(base_cols), col_name))  
})

# bin genome to 10 kb-sized "windows"
binnedgenome = bin_genome("hg38", bin_size = 10000)
binnedgenome

# matching our cnv data to the genomic bins
# necessary for plotting the individual heatmap "boxes".
# generating character matrix containing "amp", "del" or NA (i.e. no cnv in that region)

char_mat = NULL
fin = length(bc_sublists)
for (i in 1:fin) {
  bed = bc_sublists[[i]]
  gr_bed = GRanges(seqnames = bed[, 1], 
                   ranges = IRanges(
                     bed[,2], 
                     bed[,3]))
  char_mat = cbind(
    char_mat, 
    normalize_genomic_signals_to_bins(
      gr_bed, 
      bed[, 7],
      empty_value = NA))
}

# saving character matrix for making plots in second script
qsave(char_mat,
      paste0(dir, "/", ext, "/cnv_char_mat_for_heatmap.qs"),
      nthreads = 32)




