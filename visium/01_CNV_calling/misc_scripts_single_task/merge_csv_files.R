# merge multiple csv files 
library(dplyr)
library(readr)

# dir, patterns of csv files to merge
d = "/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2/whole_patient_cnvcalling/20251004/"
i = paste0(d, "01_UNMERGED") #in
o = paste0(d, "02_MERGED") # out
patterns = c("LN0025",
             "LN0027",
             "LN0438",
             "LN0193",
             "LN0040HD",
             "LN0035HD")


# define funtion
merge_csv = function(dir, pattern, output) {
  files = list.files(dir, pattern, full.names = TRUE)
  dfs = lapply(files, function(f) {
    read.table(f, header = TRUE, sep = ",")
  })
  merged = bind_rows(dfs)
  write_csv(merged, output)
}

# iter over all patient ids
lapply(patterns, function(pat) {
  out = paste0(o, "/", pat, "_merged.csv")
  merge_csv(i, pat, out)
})
