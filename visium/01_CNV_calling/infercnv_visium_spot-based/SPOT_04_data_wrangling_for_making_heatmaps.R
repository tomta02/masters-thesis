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
cnv[is.na(cnv)] = "NA" # convert NA values to string as otherwise "normalize_genomic_signals_to_bins()" function gives an error -> look into later

# PROBLEM: character matrix turned out to be way too large to compute heatmaps in R.
# replace all character values "amp", "del" by integer numbers and make mat sparse in the end.

# BIGGEER PROBLEM: if I assign integer values to "NA", "amp" and "del" now (e.g. NA = 0, amp = 1, del = 2), the 
# "normalize_genomic_signals_to_bins"function will average them to our bin size -> e.g. if I assign int num 1 to "amp", 
# and a given amplified region in a spot is e.g. 2 kb large -> 
# "normalize_genomic_signals_to_bins" with bin size 10kb will average the  2kb region of "1" (="amp") to -> 0.2 over 10kb. 
# Depending on the size of the amplified or deleted segments, 
# worried I wont be able to decipher the final matrix (what region was amp, what region was del?) since the averaged numbers 
# over our 10kb bin size might span a large range
# stupid solution: maybe assign v large integer to amp and very small to del, then differences should be large enough
# to infer the genomic state
# for now -> compute character matrix, so we can keep "amp" and "del" values, after constructing mat, convert "amp" and "del"
# to 1 and 2 (and 0 for "NA", i.e. no cnv was found) and save as sparse matrix
# think about


# the "normalize_genomic_signals_to_bins" function does not work with a full character matrix (only numeric matrices)
# -> need to supply our data as individual vectors. To do this more easily, split our cnv data per visium spot barcode
# and save thesee "slices" in list
# there's prob a way better way to do this, think about
base_cols = cnv[, 1:6] # are metadata cols, exclude

# get all the patient columns
# NOTE: extracting pat_id roi-based cnv calling - for per-patient analysis (instead of per roi), pat_id = ext
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

# old approach, did not save bc name and cbind is very slow in such a large loop - 
# try generating list first and then making df at end (see bottom approach w "charlist")

# 20251017 to finish analysis script: make small 15 bc analysis
#bc_sublists_small = bc_sublists[1:20]

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
qsave(char_mat, paste0(dir, "/", ext, "/50th_quantile/", ext, "_SPOT_04_cnv_char_mat_for_heatmap_cnv_flt_0.9.qs"), nthreads = 32)




