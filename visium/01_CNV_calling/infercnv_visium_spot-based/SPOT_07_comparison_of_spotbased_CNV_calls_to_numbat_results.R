#################################################
# SPOT_07_inferCNV_numbat_matching_spotbased
# compare cnvs found in numbat from the single-
# cell data to cnvs found per spot in the spatial
# data ("reference CNVs", script SPOT_06)
#################################################

args = commandArgs(trailingOnly = TRUE)

dir = args[1] 
nb = args[2]
infcnv = args[3] 
ext = args[4] 


# load libraries
library(qs)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(GenomicRanges)
library(VennDiagram)
library(gtrellis)
library(circlize)

# load numbat data
numbatlist = qread(nb, nthreads = 32)

# loading genotypes that were significant for numbat
gt = qread("/g/saka/Kristy/projects/composite/analysis/rdata/00_numbat_gt_forTatjana.qs")

# extract cnv calls from data
nb_cnv_list = list()

for (name in names(numbatlist)) {
  pat = numbatlist[[name]][["segs_consensus"]]
  nb_cnv_list[[name]] = pat
}

# subset nb_cnv_list by the genotypes that were significant
nb_cnv_list_subs = mapply(function(patient, genot) {
  patient[patient$seg_cons %in% genot, ]
}, nb_cnv_list, gt, SIMPLIFY = FALSE)


# now load infercnv output data
infcnv_ref = read.csv(infcnv)
infcnv_ref$chr = as.numeric(sub("chr", "", infcnv_ref$chr))

# also load information which roi
# for extracting the correct patient data from numbat data:
pat_id = sub("_.*", "", ext)

# subset numbat cnvs
nb_cnv = as.data.frame(nb_cnv_list_subs[[pat_id]])


# NOW MATCH
compare_intervals = function(df_siCNV, df_NB) {
  results = data.frame(chr = integer(),
                       state = character(),
                       match = logical(),
                       length_siCNV = integer(),
                       length_NB = integer(),
                       nr_match = integer(),
                       pct_match_siCNV = numeric(),
                       pct_match_NB = numeric())
  
  for (i in 1:nrow(df_siCNV)) {
    for (j in 1:nrow(df_NB)) {
      
      #chr matches?
      if (df_siCNV$chr[i] == df_NB$CHROM[j]) {
        
        # state matches?
        match_status = df_siCNV$state[i] == df_NB$cnv_state[j]
        
        if (match_status) {
          # make seqs based on start and end values
          seq_siCNV = seq(df_siCNV$start[i], df_siCNV$end[i])
          seq_NB = seq(df_NB$seg_start[j], df_NB$seg_end[j])
          
          # calc lengths, find intersect, calc percentages
          length_siCNV = length(seq_siCNV)
          length_NB = length(seq_NB)
          
          seq_intersect = intersect(seq_siCNV, seq_NB)
          nr_match = length(seq_intersect)
          if (length(seq_intersect) > 0) {
            start_match = min(seq_intersect)
            end_match = max(seq_intersect)
          } else {
            start_match = NA
            end_match = NA
          }
          
          pct_match_siCNV = nr_match / length_siCNV * 100
          pct_match_NB = nr_match / length_NB * 100
          } else {
            start_match = NA
            end_match = NA
            length_siCNV = NA
            length_NB = NA
            nr_match = NA
            pct_match_siCNV = NA
            pct_match_NB = NA
            }
        
        # apend to results
        results = rbind(results, data.frame(chr = df_siCNV$chr[i],
                                            state = df_siCNV$state[i],
                                            match = match_status,
                                            length_siCNV = length_siCNV,
                                            length_NB = length_NB,
                                            nr_match = nr_match,
                                            start_match = start_match,
                                            end_match = end_match,
                                            pct_match_siCNV = pct_match_siCNV,
                                            pct_match_NB = pct_match_NB))
      }
    }
  }
  
  return(results)
}

shared = compare_intervals(infcnv_ref, nb_cnv)
shared = shared[shared$match == TRUE & shared$nr_match > 0,]

# make plots for reference cnvs, numbat cnvs, and shared
# first step: converting data.frames to granges objects
make_granges = function(df, chr_col, start_col, end_col, state_col) {
  GRanges(
    seqnames = df[[chr_col]],
    ranges = IRanges(start = df[[start_col]], end = df[[end_col]]),
    cnv_state = df[[state_col]]
  )
}


gr_infcnv = make_granges(infcnv_ref, "chr", "start", "end", "state")
gr_numbat = make_granges(nb_cnv, "CHROM", "seg_start", "seg_end", "cnv_state_post")
gr_shared = make_granges(shared, "chr", "start_match", "end_match", "state")

seqlevelsStyle(gr_infcnv) = "UCSC"
seqlevelsStyle(gr_numbat) = "UCSC"
seqlevelsStyle(gr_shared) = "UCSC"


pdf(
  paste0(dir, "/", ext, "/SPOT_07_infercnv_numbat_comparison.pdf"),
  width = 22,
  height = 6
  )


col_fun = c("amp" = "red", "del" = "blue")


gtrellis_layout(
  category = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"),
  species = "hg38",
  add_name_track = TRUE,
  name_fontsize = 12,
  add_ideogram_track = TRUE,
  n_track = 3, 
  equal_width = FALSE,
  title = paste("CNVs inferred using inferCNV, numbat, and overlapping CNVs - ", ext),
  title_fontsize = 22,
  track_ylab = c("inferCNV", "numbat", "shared"),
  lab_fontsize = 18,
  track_axis = FALSE,
  xlab = "")

add_rect_track(
  gr_infcnv,
  track = 2,
  h1 = 1, 
  h2 = 0,
  gp = gpar(col = NA, fill = col_fun[mcols(gr_infcnv)$cnv_state])
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
