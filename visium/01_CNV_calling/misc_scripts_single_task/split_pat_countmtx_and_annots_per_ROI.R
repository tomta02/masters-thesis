# running inferCNV on all visium datasets grouped per patient in a "per-spot" 
# resolution is taking way too long (bigmem queue, max mem available, 32 threads, 
# multiple days & counting) -> splitting up count matrix and annot file per
# ROI (i.e. individual visium slides, not grouped per patient) instead.

args = commandArgs(trailingOnly = TRUE)

vislist_path = args[1]
annot_path = args[2]
countmtx_path = args[3]
outdir_path = args[4]

# load libraries
library(qs)
library(tidyverse)

# load list of preprocessed visium seurat objects, extract their names (i.e names of all ROI)
visiumlist = qread(vislist_path, nthreads = 32)
roinames = names(visiumlist)
rm(visiumlist)

# load patient annot files & count matrices
# they will be subset based on ROI names from "roinames" list
annot = read.table(
  file = annot_path,
  sep = "\t",
  header = FALSE,
  col.names = c("barcode", "spatial_domain")
  )

countmtx = read.table(
  file = countmtx_path,
  sep = "\t",
  header = TRUE
)

countmtx = column_to_rownames(countmtx, "Genes") %>% t


# for subsetting count mtx and annot: need to loop through all roi names
for (roi in roinames) {
  # row names start with 'roi', followed by underscore
  match_rows_cntmtx <- grep(paste0("^", roi, "_"), rownames(countmtx), value = TRUE)
  # in annot df, rois are saved in col "barcode"
  match_rows_annot_idx <- grep(paste0("^", roi, "_"), annot$barcode)
  
  # we also need to keep the barcodes of the control group for both the
  # count mtx and the annot file... they don't have underscores:
  no_underscore_bcs_countmat = grep("_", rownames(countmtx), invert = TRUE, value = TRUE)
  no_underscore_bcs_annot_idx = grep("_", annot$barcode, invert = TRUE)
  
  # concatenate 
  all_bc_to_keep_cntmat = c(match_rows_cntmtx, no_underscore_bcs_countmat)
  all_bc_to_keep_annot = c(match_rows_annot_idx, no_underscore_bcs_annot_idx)
  
  if (length(match_rows_cntmtx) > 0) {
    # subset count mat
    subset_df <- countmtx[all_bc_to_keep_cntmat, , drop = FALSE] %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("Genes")
    
    # only in the countmatrix the barcodes have a "-" in the name, in the annot they have a "."
    # replace the "." in the countmat barcode by "-"
    colnames(subset_df)[-1] = gsub("\\.", "-", colnames(subset_df)[-1])
    
    # subset annot df by row indices extracted earlier
    subset_annot = annot[all_bc_to_keep_annot, , drop=FALSE]
    
    # asign to new var in the global environment, named after roi
    assign(roi, subset_df)
    assign((paste0(roi, "_annot")), subset_annot)
    
    write.table(
      subset_df,
      file = paste0(outdir_path, roi, "_countmtx.tsv"),
      sep = "\t",
      row.names = FALSE,
      #quote = FALSE
    )
    
    write.table(
      subset_annot,
      file = paste0(outdir_path, roi, "_annot.tsv"),
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
    
  } else {
    # print in case no matches in count mtx
    message("NOTE: no matching entries found for ROI: ", roi)
  }
}



