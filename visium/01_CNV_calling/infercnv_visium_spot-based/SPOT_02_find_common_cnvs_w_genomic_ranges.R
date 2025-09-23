#########################################
# Cleaning CNV data with GenomicRanges #
#########################################

args = commandArgs(trailingOnly = TRUE)

dir = args[1]
cnvs = args[2]
ext = args[3]


# the CNVs identified per individual Visium spot are quite noisy, and we need 
# to identify the most confident CNVs. Done using the GenomicRanges library 
# load libraries
library(tidyverse)
library(magrittr)
library(qs)
library(GenomicRanges)
library(tidyverse)

# load sample
cnvs_to_reformat = read.table(cnvs, header = TRUE)

# subset cnvs by cnvs found in patient, not cnvs found in negative control 
pat_id = sub("_.*", "", ext)
bc_to_keep = grep(paste0("^", pat_id, "_"), cnvs_to_reformat$bc_name)
cnv_flt = cnvs_to_reformat[bc_to_keep, , drop=FALSE]

gr = GRanges(seqnames = cnv_flt$chr,
             ranges = IRanges(
               start = cnv_flt$start,
               end = cnv_flt$end,
               names = cnv_flt$cnv_name),
             barcode = cnv_flt$bc_name,
             state = cnv_flt$state_mcmc,
             cnv_name = cnv_flt$cnv_name
)

# disjoin genomic ranges -> divide sequences found into non-overlapping "chunks" 
# count how often each of the non-overlapping segments appears, 
# take 75% quantile (i.e. most common segments across most barcodes) to keep high
# confidence CNVs and get rid of noise

# disjoin. Reverse mapping = true, to keep metadata
dj = disjoin(gr, with.revmap = TRUE)

mcols(dj) = do.call(rbind, 
                    lapply(mcols(dj)$revmap, function(i){
                      data.frame(
                        bc = paste(mcols(gr)$barcode[i], collapse=","),
                        state = paste(mcols(gr)$state[i], collapse=",")
                      )
                    })
)

# count how often each sequence chunk appears in "gr", by using "countOverlaps"
ol = countOverlaps(dj, gr)
mcols(dj)$count = ol

# clean up disjoined genomic ranges: remove chunks that appear with very low frequency (i.e. noise) across all chromosomes
# calc 75th percentile (i.e. top 25%) threshold
thresh = quantile(mcols(dj)$count, 0.75, na.rm = TRUE)
dj_clean = dj[mcols(dj)$count >= thresh]


# neext, expand cleaned disjoined fragments (dj_clean), to have each barcode in separate row
# extract metadata
expanded_list = lapply(seq_along(dj_clean), function(i) {
  bc_vec = unlist(strsplit(mcols(dj_clean[i])$bc, ","))
  state_vec = unlist(strsplit(mcols(dj_clean[i])$state, ","))
  
  # repeat granges row [i] for each bc
  n = length(bc_vec)
  gr = dj_clean[i]
  gr_expanded = rep(gr, n)
  
  # update netadata
  mcols(gr_expanded)$bc = bc_vec
  mcols(gr_expanded)$state = state_vec
  
  return(gr_expanded)
})

# combine all into 1 GRanges
dj_clean_expanded = do.call(c, expanded_list)

# now, convert expanded granges to bed-style format - needed for most plotting libraries
# first, convert it to df
gr_df = data.frame(
  "chr" = dj_clean_expanded@seqnames,
  "start" = dj_clean_expanded@ranges@start,
  "end" = (dj_clean_expanded@ranges@start + dj_clean_expanded@ranges@width),
  "width" = dj_clean_expanded@ranges@width,
  "width_kb" = (dj_clean_expanded@ranges@width)/1000,
  "strand" = dj_clean_expanded@strand,
  "barcode" = dj_clean_expanded$bc,
  "state" = dj_clean_expanded$state
)

# now: reshape gr_df into sparse BED-like format, needed for most plotting libraries
# ISSUE: Issue with disjoining and regrouping genomic regions: 
# across different visium spots, the "same" copy number variation (i.e. same genomic region, 
# same underlying CNV: "amp", "del") might have a slight difference in state: 
# e.g. the same genomic region might be predicted to be amplified or deleted across many spots, 
# but in some spots the region might have been identified as state "4" (i.e. amplification of singular allele),
# whereas in other spots it might have been called as state "5" (amplification of both alleles). 
# -> when cleaning up the cnv df to be exported, this creates issues - instead of singular 
# integer state (e.g. "4") we used to report multiple states concatenated as string, separated by comma (e.g. "4,5"). 
# Resolve this with a majority vote below, to only report singular states

# write function to resolve multiple cnv states
# very clumsy but only "adjacent" dual states were reported: no "triple states" (e.g. "4,5,6") and no unrelated states (e.g. "1,4")
# therefore save these as "state_pairs_ok" and raise a warning in case the reported cnv state does not match "state_pairs_ok"
stateresolve = function(x, 
                        state_pairs_ok = c("1,2", "2,1", # are both deletions
                                           "4,5", "5,4", # are both amplifications
                                           "5,6", "6,5")) { # are both amplifications 
  x = na.omit(x) 
  
  # if state = singular value: return as is, no problem
  if (length(x) == 1) return(x)
  
  # if we have multiple vals like ("4", "4, "4", "5") then keep only unique ones ("4", "5") and collapse them
  s = paste(unique(x), collapse = ",")
  
  # if state matches an "allowed" paiir, do max vote
  if (s %in% state_pairs_ok) {
    vals = unlist(strsplit(s, ","))
    winner = names(sort(table(vals), decreasing = TRUE))[1]
    return(winner)
  }
  
  # in case unexpected state combi occurs, warn and return NA
  warning(glue::glue("weird state combination: {s}"))
  return(NA_character_)
}


# now: reshape gr_df into sparse BED-like format, needed for most plotting libraries - use "stateresolve" function to take care of conflicting state issues
gr_df_wide <- gr_df %>%
  pivot_wider(
    names_from = barcode,
    values_from = state,
    values_fn = list(state = ~ stateresolve(.)),
    values_fill = NA
  )

# save finalized object
write.csv(gr_df_wide,
file = paste0(dir, "/", ext, "/SPOT_02_cnvs_cleaned.csv"),
row.names = FALSE)

# create 2nd df with binarized CNV state: "amp" or "del" instead of CNV states 1-6
# (see: https://github.com/broadinstitute/infercnv/wiki/infercnv-i6-HMM-type for more info on cnv states)

gr_df_wide_bin = gr_df_wide %>%
  mutate(across(7:ncol(.), ~ case_when(
    . %in% c("1", "2") ~ "del",
    . %in% c("4", "5", "6") ~ "amp",
    TRUE ~ as.character(.)
  )))


write.csv(gr_df_wide_bin,
          file = paste0(dir, "/", ext, "/SPOT_02_cnvs_cleaned_binarized.csv"),
          row.names = FALSE)
