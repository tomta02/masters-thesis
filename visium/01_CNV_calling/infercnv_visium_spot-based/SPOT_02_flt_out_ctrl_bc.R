#########################################
# Cleaning CNV data with GenomicRanges #
#########################################

args = commandArgs(trailingOnly = TRUE)

dir = args[1]
cnvs = args[2]
ext = args[3]

# subset cnvs by cnvs found in patient, not cnvs found in negative control - analyzed separately, is just noise

library(tidyverse)
library(magrittr)
library(qs)
library(GenomicRanges)
library(tidyverse)

# load sample
cnvs_to_reformat = read.table(cnvs, header = TRUE)

# get patient IDs
pat_id = sub("_.*", "", ext)
bc_to_keep = grep(paste0("^", pat_id, "_"), cnvs_to_reformat$bc_name)
cnv_flt = cnvs_to_reformat[bc_to_keep, , drop=FALSE]

write.csv(cnv_flt,
           file = paste0(dir, "/", ext, "/NEW_SPOT_02_cnv_prob_greaterthan_0.9_bc_flt_no_normal_bc.csv"),
           row.names = FALSE)

