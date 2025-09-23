args = commandArgs(trailingOnly = TRUE)

cnv_regs = args[1]
mcmc_obj = args[2]
out_dir = args[3]
ext = args[4]

# loading libraries
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(magrittr)
library(qs)
library(infercnv)

set.seed(99)

# create output directory
outs_path = paste0(out_dir, ext)
dir.create(outs_path)


# loading output of "spot-based" infercnv run, singular ROI: LN0438_MAAFHY1_R1
# raw results (post MCMC)
print("[x] reading in CNV regions table.")
cnv_reg = read.table(cnv_regs, sep = "\t", header = TRUE)

# read in whole Visium spot-based mcmc infercnv object
print("[x] reading in mcmc object.")
mcmc = readRDS(mcmc_obj)


# start of data reformating: "cell_probabilities" slot in mcmc object.
print("[x] extracting necessary properties from mcmc object.")
cellprobs = mcmc@cell_probabilities

# get their names from cnv_regions slot
cnvregs = mcmc@cnv_regions
names(cellprobs) = as.character(cnvregs)

# cellprobs not state the barcode where each cnv appears - get bc information from count.data, and match it to cnvs in second step using the integer barcode indices in "cell_gene" slot of mcmc object
bc = colnames(mcmc@count.data)

# obtain barcode names by matching bc indices to bc names
cellprobs = map2(
  cellprobs,
  mcmc@cell_gene,
  ~ {
    # integer barcode indices from cell_gene, look for bc names in "bc"
    barcode_idx = as.integer(.y$Cells)
    colnames(.x) = bc[barcode_idx]
    .x
  }
)

# use "cellprobs" list to make df of cnvs - to later match match with raw output file "cnv_regs_spots_r" to get remaining cols that are necessary for downstream analysis (chr, start & stop genomic coords)
print("[x] reformat CNV information.")

cnvs_raw = data.frame(
  cnv_name = names(cellprobs),
  bc_name = sapply(cellprobs, function(x) colnames(x)),
  p = sapply(cellprobs, function(x) max(x[, 1])),
  state_mcmc = sapply(cellprobs, function(x) rownames(x)[which.max(x[, 1])]),
  stringsAsFactors = FALSE
)

# prototype of final res: extracted info from mcmc object (= post mcmc, post bayes net -> only keep most probable cnvs) & merged with most important cols of raw "cnv_regs_spots_r" output file
cnvs_proto = cnvs_raw %>%
  left_join(
    cnv_reg %>%
      select(cnv_name, cell_group_name, state, chr, start, end),
    by = c("cnv_name" = "cnv_name", "bc_name" = "cell_group_name")
  )

# last steps: apply qc filtering thresholds to cnvs called: remove all cnvs with reported state "3" in col "state_mcmc" (="normal", no amp/del), remove all cnvs with probability (col: "p") < 0.7, to keep only high confidence cnvs

cnvs_fin = cnvs_proto[cnvs_proto$state_mcmc != 3, ] #&
#cnvs_proto$p > 0.7, ]
# changed approach: export all cnvs, no matter the probability, and do filtering based on CNV probability later downstream.

# additionally, add some further columns to "cnvs_fin" that will be useful downstream
cnvs_fin$width = cnvs_fin$end - cnvs_fin$start
cnvs_fin$width_kb = (cnvs_fin$width)/1000


# save final result
print("[x] saving final results.")
write.table(cnvs_fin, 
            paste0(outs_path, "/cnv_final.csv"),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

# save raw data as well - since whole mcmc object is too large to save it whole, and is currently saved on scratch (no backups)
# save the components used to create the "cnv_final.csv" in case it ever needs to be built again
write.table(cellprobs, 
            paste0(outs_path, "/RAW_cell_probabilities_cnv.csv"))

write.table(cnvregs,
            paste0(outs_path, "/RAW_cnv_regions.csv"),
            row.names = TRUE,
            col.names = FALSE)

write.table(bc,
            paste0(outs_path, "/RAW_spot_barcode_names.csv"),
            row.names = TRUE,
            col.names = FALSE)
