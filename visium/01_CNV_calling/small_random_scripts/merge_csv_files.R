# merge multiple csv files 
library(dplyr)
library(readr)

# dir, patterns of csv files to merge
d = "/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2/20251004/01_CNVs_GROUPED_PER_PATIENT/"
patterns = c("LN0025",
             "LN0027",
             "LN0438",
             "LN0193",
             "LN0040HD",
             "LN0035HD")


# define funtion
merge_csv = function(dir, pattern, output) {
  files = list.files(dir, pattern, full.names = TRUE)
  dfs = lapply(files, read_csv)
  merged = bind_rows(dfs)
  write_csv(merged, output)
}

# iter over all patient ids
lapply(patterns, function(pat) {
  out = paste0(d, "/", pat, "_merged.csv")
  merge_csv(d, pat, out)
})



