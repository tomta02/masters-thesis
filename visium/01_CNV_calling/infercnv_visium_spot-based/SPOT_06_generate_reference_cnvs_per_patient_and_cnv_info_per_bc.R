############################################################################
# SPOT_06: generating reference set of cnvs per patient and identifying
# presence/absence of each reference cnv per each barcode via state maxvote.
# while this is step 6 of "CNV_calling", the outputs of this script are 
# necessary for for correlating pathway activity per spot with the presence
# or absence of certain CNVs (see "03_CNVxPW_analysis")
############################################################################

args = commandArgs(trailingOnly = TRUE)

dir = args[1] #"/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2/roi_based_cnvcalling"
cnvs = args[2] #paste0(dir, "/", ext, "/SPOT_02_cnvs_cleaned_binarized.csv") 
ext = args[3] #"LN0438_MAAFHY1_R1" 
outdir = args[4] #paste0("/g/saka/Tatjana/data/03_CNVxPW_analysis/", ext) 


# load libraries
library(tidyverse)
library(magrittr)
library(qs)
library(GenomicRanges)
library(tidyverse)
library(gtrellis)

# load cnv data for LN0438_MAAFHY1_R1
cnv = read.csv(cnvs)
cnv[is.na(cnv)] = "NA"

# we are reducing cnvs across all barcodes to generate a "reference" set of cnvs
# see https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html
# takes set of overlapping genomic intervals and merges them into the smallest set of non-overlapping 
# intervals that cover same genomic positions. 
# e.g. in: [1,5], [3,7], [10,12] -> out: [1,7], [10,12] (we find union of cnvs across all barcodes, to identify
# consensus cnvs, since otherwise the data is too noisy)
# for this, need to convert cnvs to "GenomicRanges" object - first "pivot_longer", then convert to gr
cnv_long = cnv %>%
  pivot_longer(
    cols = 7:ncol(.),       
    names_to = "barcode",
    values_to = "state"    
  )

gr = GRanges(seqnames = cnv_long$chr,
             ranges = IRanges(
               start = cnv_long$start,
               end = cnv_long$end),
             barcode = cnv_long$barcode,
             state = cnv_long$state
)

# MAKING 'REFERENCE' cnvs
ref = reduce(gr, with.revmap = TRUE)

mcols(ref) = do.call(rbind, 
                     lapply(mcols(ref)$revmap, function(i){
                       data.frame(
                         barcode = paste(mcols(gr)$barcode[i], collapse=","),
                         state = paste(mcols(gr)$state[i], collapse=",")
                       )
                     }))

# next step: now that genomic ranges object was reduced, we need to decide whether a given genomic range
# counts as "amplified" or "deleted" - there will be a clear trend (i.e. if a cnv is actually there, it will have been
# called as "amp" or "del" in the majority of barcodes - but data is highly noisy and for some bc the opposite is predicted
# i.e. "del" when the overarching trend is "amp" or the other way around)
# for this -> do majority vote

# to get state of reference cnvs:
get_state_ref = function(x) {
  # states from reverse mapped metadata of "ref"
  states = strsplit(x, ",", fixed = TRUE)[[1]]
  states = states[states %in% c("amp", "del")]
  
  # count how many amp or del
  n_amp = sum(states == "amp")
  n_del = sum(states == "del")
  
  # max vote
  if (n_amp > n_del) {
    return("amp")
  } else if (n_del > n_amp) {
    return("del")
  } else {
    return("tie")
  }
}

# for making majority vote, we need to convert reference cnvs from genomic ranges object
# back into a data.frame. 
# first, get majority vote to decide the states of the reference cnvs, then add this as col to the new df

ref_states = lapply(mcols(ref)$state, get_state_ref)

ref_df = data.frame(
  chr = as.character(seqnames(ref)),
  start = start(ref),
  end = end(ref),
  width = width(ref),
  width_kb = width(ref) / 1000,
  state = unlist(ref_states)
)

# T_T to make neat plots of the cnvs found in each patient, we need to convert the data back into a genomicranges
# object, since the gtrellis library works with data GenomicRanges format... 

ref_gr_fin = GRanges(seqnames = ref_df$chr,
                     ranges = IRanges(
                       start = ref_df$start,
                       end = ref_df$end),
                     state = ref_df$state
)


# making gtrellis plot
out_pdf = paste0(outdir, "/", ext, "_reference_cnv_states.pdf")
pdf(out_pdf, width = 20, height = 2)

col_fun = c("amp" = "red", "del" = "blue")

gtrellis_layout(
  species = "hg38",
  add_name_track = TRUE,
  add_ideogram_track = TRUE,
  n_track = 1,
  equal_width = FALSE,
  title = paste0("reference CNVs across all Visium spots - ", ext),
  #track_ylab = c(""),
  track_axis = FALSE,
  xlab = "")

add_rect_track(
  ref_gr_fin,
  track = 2,
  h1 = 1,
  h2 = 0,
  gp = gpar(col = NA, fill = col_fun[mcols(ref_gr_fin)$state])
)
dev.off()


# NEXT STEP
# maybe split script in 2 parts here?
# in order to determine for each individual barcode whether it has a given cnv or not, 
# we need to apply GenomicRanges::reduce to all intra-barcode chunks 
# this is due to how the data was processed & cleaned in script "SPOT_03_find_common_cnvs_w_genomic_ranges.R"
# in order to get rid of noise.
# in essence, for each spot, we have exactly the same nr. and sizes of genomic bins present. 
# most of these will have been called with the same cnv_state ("amp", "del") like the cnvs of our reference dataset
# but some not. In order to figure out which spots harbour the same cnvs we have in our reference cnvs:
# 1.) convert cnvs into GenomicRanges object and apply intra-barcode GenomicRanges::reduce (with collecting "state" info by revmap)
# 2.) onduct majority vote: is a given cnv present or absent?, convert data back to data.frame
# 3.) Compare the "state" of the cnvs per barcode identified in step "2.)" with reference cnvs and perform matching operation:
#           -> if cnv in barcode did not pass max vote (more ranges with "NA" compared to "amp"/"del") -> cnv absent, cnv_test = 0
#           -> if cnv in barcode matches reference cnv state for that patient -> cnv_test = 1
#           -> if cnv in barcode does not match reference cnv state for that patient -> cnv_test = 2

# 1.) converting cnvs into genomicranges object (per barcode)"
# first, pivot_longer to get data in right shape

cnv_long_df = data.frame(cnv_long)
cnv_long_df = na.omit(cnv_long_df)
cnv_spots = split(cnv_long_df, cnv_long_df$barcode)

# make genomic ranges list (1 barcode = 1 list element)
gr_lst = lapply(cnv_spots, function(df) {
  GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start, end = df$end),
    barcode = df$barcode,
    state = df$state
  )
})


# apply GenomicRanges::reduce whiel keepin track of state data
red_gr_lst = lapply(gr_lst, function(gr_obj) {
  
  red = reduce(gr_obj, with.revmap = TRUE)
  
  # keep metadata using revmap., saving state and bc info
  aggregated_mcols = do.call(rbind, lapply(mcols(red)$revmap, function(i) {
    data.frame(
      barcode = paste(mcols(gr_obj)$barcode[i], collapse = ","),
      state = paste(mcols(gr_obj)$state[i], collapse = ",")
    )
  }))
  mcols(red) = aggregated_mcols
  
  red
})

# 2.) like above: for each barcode (i.e. each element of "red_gr_list"), do a majority voting to get CNV state:
# make func for this. This function slightly different from the one "get_ref_bc", bc it includes the possibility
# of a cnv being absent "NA" in certain barcodes, this is obv not possible for the reference cnvs
# revisit to make singular function appliccable to both cases later
get_state_bc = function(x) {
  states = strsplit(x, ",", fixed = TRUE)[[1]]
  states = states[states %in% c("amp", "del", "NA")]
  
  n_amp = sum(states == "amp")
  n_del = sum(states == "del")
  n_na  = sum(states == "NA")
  
  if (n_na > n_amp && n_na > n_del) {
    return(NA)
  } else if (n_amp > n_del) {
    return("amp")
  } else if (n_del > n_amp) {
    return("del")
  } else {
    return("tie")
  }
}


bc_states = lapply(red_gr_lst, function(gr) {
  lapply(mcols(gr)$state, get_state_bc)
}
)

# now, convert GenomicRanges object of each barcode into data.frame, save in list of data.frames
gr_df_list = lapply(seq_along(red_gr_lst), function(i) {
  
  gr = red_gr_lst[[i]]
  states = bc_states[[i]]
  
  data.frame(
    chr = as.character(seqnames(gr)),
    start = gr@ranges@start,
    end = (gr@ranges@start + gr@ranges@width),
    width = gr@ranges@width,
    width_kb = width(gr) / 1000,
    strand = as.character(strand(gr)),
    state = unlist(states) 
  )
})

names(gr_df_list) = names(red_gr_lst)

# across all barcodes, the reduced ranges identified should be exactly the same as in the reference
# object: just not clear if they passed the max vote and actually count as "cnv"
# double check whether ranges are the same across all bc

cols_to_check = c("chr", "start", "end")

colslist = lapply(gr_df_list, function(df) df[cols_to_check])
allsame = all(sapply(colslist, function(x) identical(x, colslist[[1]])))

allsame
# TRUE 

# 3.) matching operation: 
# if gr_df_list[[barcode]]$state = NA (i.e. CNV is absent), cnv_test = 0
# if gr_df_list[[barcode]]$state matches ref_df$state, cnv_test = 1
# else (i.e mismatch between gr_df_list[[barcode]]$state and ref_df$state), cnv_test = 2

gr_df_list = imap(gr_df_list, ~ .x %>%
                    mutate(
                      cnv_test = case_when(
                        is.na(state) ~ 0L,
                        state == ref_df$state ~ 1L,
                        TRUE ~ 2L),
                      # fill col "bc" with name(gr_df_list[bc]), to keep track of barcode name
                      bc = .y
                    ))

# each df in"gr_df_list" correspond to a different barcode, and within each df/each barcode, each row corresponds to 
# a different cnv found within that barcode
# perform row-wise stacking: right now we "per barcode view": singular barcode and all cnvs for that barcode; 
# but we want "per cnv view": singular cnv and all barcodes + info whether a given bc has a cnv or not. 

colstokeep = c("chr", "start", "end", "cnv_test", "bc")

# num of rows for each df should be exactly the same, check again
# make vector containing nr. of rows of each df as elements
nrows_v = map_int(gr_df_list, nrow)

if (!all(nrows_v == nrows_v[1])) {
  stop("dfs have diff nr of rows",
       paste(nrows_v, collapse = ", "))
}

n_rows = nrows_v[1]

# list with each element corresponding to a given row index (rows = the cnvs)
cnv_list = map(1:n_rows, function(i) {
  # for the current row index: extract *that specific* row from all dfs 
  map_dfr(gr_df_list, function(df) {
    df[i, colstokeep]
  })
})


# add unique index to cnv, to keep track of cnvs later 
cnv_list = Map(function(df, i) {
  df$cnv_id = i #paste0("cnv_", i)
  df
}, cnv_list, seq_along(cnv_list))


## RESULT: cnv_list contains following info: 
# 1.) length of list corresponds to nr. of cnvs that was found,
# 2.) for each cnv we have information on which barcodes the cnv was present in,
# and which barcodes the cnv was not found in. (based on majority voting of genomic "chunks"
# that were obtained from cleaning the data in script "SPOT_03_find_common_cnvs_w_genomic_ranges.R")

# save this list for next step - comparing pathway activities found per barcode with the presence or absence of certain cnvs
qsave(cnv_list, paste0(outdir, "/", ext, "_01_list_of_all_cnvs_found.qs"))
