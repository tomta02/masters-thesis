#########################################
# data wrangling for heatmap generation #
#########################################

args = commandArgs(trailingOnly = TRUE)

dir = args[1] # "/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2/roi_based_cnvcalling" #
cnvs = args[2] # paste0(dir, "/", ext2, "/NEW_SPOT_06_cleaned_cnvs_per_bc_from_NEW_SPOT_03_cnvs_cleaned_binarized_for_mking_charmat.csv") #
ext = args[3] # "LN0438_MAAFHY1_R1" #



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
cnv[is.na(cnv)] = "NA" # convert NA values to string as otherwise "normalize_genomic_signals_to_bins()" function gives an error -> look into later

# the "normalize_genomic_signals_to_bins" function does not work with a full character matrix (only numeric matrices)
# -> need to supply our data as individual vectors. To do this more easily, split our cnv data per visium spot barcode
# and save thesee "slices" in list
# there's prob a way better way to do this, think about
base_cols = cnv[, 1:6] # are metadata cols, exclude

# get all the patient columns
# NOTE: extracting pat_id for roi-based cnv calling (not appliccable when analyzing just ROIs)
#pat_id = sub("_.*", "", ext)

bc_cols = grep(ext, names(cnv), value = TRUE)

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

# old approach did not save bc name and cbind is very slow in such a large loop - 
# try generating list first and then making df at end (see bottom approach w "charlist")

fin = length(bc_sublists)
charlist = vector("list", fin)
colnames = character(fin)

for (i in seq_len(fin)) {
  bed = bc_sublists[[i]]
  gr_bed = GRanges(
    seqnames = bed[, 1],
    ranges = IRanges(bed[, 2], bed[, 3])
  )
  mat_i = normalize_genomic_signals_to_bins(
    gr_bed,
    bed[, 7],
    empty_value = NA
  )
  charlist[[i]] = mat_i
  colnames[i] = names(bed)[7]
}

char_mat = do.call(cbind, charlist)
colnames(char_mat) = colnames



# saving character matrix for making plots in second script
qsave(char_mat, paste0(dir, "/", ext, "/NEW_SPOT_04_char_mat_from_NEW_SPOT_03_75thqn_bin_flt_0.9.qs"), nthreads = 32)


