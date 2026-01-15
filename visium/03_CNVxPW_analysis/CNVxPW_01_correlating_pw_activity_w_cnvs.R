###########################################################################
# CNV x PW analysis
# trying to infer whether the activity of certain biological pathways is 
# associated with the presence of certain cnvs.
###########################################################################

args = commandArgs(trailingOnly = TRUE)

pw_act = args[1] # "/g/saka/Tatjana/data/02_pathway_analysis/progeny_decoupler/sct_normalized/pr_custom/pw_act_vals_and_pvals/01_pw_act_vals_all_roi_flt_pval_smaller_0.05.qs"
cnvs_new = args[2] # CNVS NEW "/g/saka/Tatjana/data/03_CNVxPW_analysis/LN0438_MAAFHY1_R1/LN0438_MAAFHY1_R1_01_list_of_all_cnvs_found_NEW_shouldbesameasold.qs"
ext = args[3] # "LN0438_MAAFHY1_R1" # 
outdir = args[4] #paste0("/g/saka/Tatjana/data/03_CNVxPW_analysis/", ext) # 

# load libraries
library(tidyverse)
library(magrittr)
library(qs)
library(tidyverse)
library(gtrellis)
library(rstatix)

cnv_list = qread(cnvs_new, nthreads = 32)

# act_flt = pathway activity values already filtered by p-value, to keep activity value only for spots that had sign. pw association
act_flt = qread(pw_act,
                nthreads = 32)

# get patient id for subsetting "act_flt" list (it contains pathway activity )
pat_id = sub("_.*", "", ext)

# subset list by patient, and then filter to only keep the barcodes belonging to current ROI
pat_act = act_flt[[pat_id]]
ext_act = pat_act[grepl(ext, rownames(pat_act)), ]

# for some pathways, there was only a marginal contribution of that pw to expression profile of certain visium spots
# meaning that there are only few visium spots associated with that pathway - disregard pathways that were only active in
# less than 10% of all visium spots. 
# For this, calc 10% threshold
th = ceiling(nrow(ext_act) * 0.1) 

# filter out pathways = columns according to scheme above
ext_act_flt = ext_act[, colSums(!is.na(ext_act)) >= th]

# replace dot by "-" in our barcode names, as below we do matching by bc - and the barcodes in the act df have a "-" instead of "."
cnv_list = map(cnv_list, ~ .x %>%
                 mutate(bc = gsub("\\.", "-", bc))
)

# iterate over cnv_list and join with LN_act_df by bc
# for this, create col "bc" in our pathway activity levels data.frame
ext_act_flt = ext_act_flt %>%
  rownames_to_column(var = "bc")


cnv_list = map(cnv_list, ~ .x %>%
                 left_join(ext_act_flt, by = "bc")
)

# stats: run shapiro test for normality for each of the pathway activity lvls (looking at each of the "cnv_test" groups individually: 
# 0, 1, 2), to determine whether to follow up w a regular unpaired 2-sample t test or a Mann Whitney U test. 
# running test on singular cnv, as most likely dist will not be normal anyway and will have to go forth with Mann Whitney U
test_cnv = cnv_list[[1]]

# define columns of interest
pw = c("Androgen",
       "Estrogen",
       "EGFR",
       "Hypoxia",
       "JAK-STAT",
       "MAPK",
       "NFkB",
       "p53",
       "PI3K",
       "TGFb",
       "TNFa",
       "Trail",
       "VEGF",
       "WNT")

# intersect our columns w pathways given in pw
cols_to_test = intersect(pw, colnames(test_cnv))

# normtest = test_cnv %>%
#   pivot_longer(
#     cols = all_of(cols_to_test),
#     names_to = "variable",
#     values_to = "value"
#   ) %>%
#   group_by(cnv_test, variable) %>%
#   summarise(
#     shapiro = list(shapiro.test(value)),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     statistic = map_dbl(shapiro, "statistic"),
#     p_value   = map_dbl(shapiro, "p.value")
#   ) %>%
#   select(-shapiro)

# since assumption of normality not met, calculate Mann Whitney U-test (Wilcox rank sum)
wilcox_list = list()

for (i in seq_along(cnv_list)) {
  df = cnv_list[[i]]
  
  # keep unique chr, start,end and cnv id
  chr_val = unique(df$chr)
  start_val = unique(df$start)
  end_val = unique(df$end)
  cnv_val = unique(df$cnv_id)
  
  wilcox_list[[i]] = df %>%
    filter(cnv_test %in% c(0,1)) %>%
    pivot_longer(
      cols = all_of(cols_to_test),
      names_to = "variable",
      values_to = "value"
    ) %>%
    group_by(variable) %>%
    wilcox_test(value ~ cnv_test) %>%
    left_join(
      df %>%
        filter(cnv_test %in% c(0,1)) %>%
        pivot_longer(
          cols = all_of(cols_to_test),
          names_to = "variable",
          values_to = "value"
        ) %>%
        group_by(variable) %>%
        wilcox_effsize(value ~ cnv_test, ci = FALSE),
      by = c("variable", "group1", "group2")
    ) %>%
    ungroup() %>%
    mutate(
      chr = chr_val,
      start = start_val,
      end = end_val,
      cnv = cnv_val
    )
}


# for plotting:
# make 1 df out of list of dfs and log-transform p val
all_wilcox_df = bind_rows(wilcox_list) %>%
  mutate(
    logp = -log10(p),
    cnv_name = paste0(chr, "_", round(start/1000), "_", round(end/1000)),
    sig = ifelse(p < 0.05, "significant", "ns"),
    variable = str_remove(variable, "\\.x$"),
    # problem: chromosomes in cnv_names are sorted alphanumerically in plot labels
    # e.g chr1, chr10, chr11, etc; instead of chr1, chr2,c hr3, ...
    # extract numeric part of each chromosome, order them, then convert to factor for them to stay that way
    chr_num = as.numeric(str_remove(chr, "chr"))  
  ) %>%
  arrange(chr_num, start) %>%
  mutate(
    cnv_name = factor(cnv_name, levels = unique(cnv_name))
  ) %>%
  select(-chr_num)  


# plot heatmap: plotting effect size in color and p value using alpha value (color
# saturated when significant, see-through when not significant
# maybe think about more eleegant way of plotting)

# also add additional identifyers of heatmap: dot if effsize > 0.3, 
# and a star if effsize > 0.5

p = ggplot(all_wilcox_df, aes(x = cnv_name, y = variable, fill = effsize, alpha = sig)) +
  geom_tile(color = "grey90", size = 0.1) +
  #scale_fill_viridis(option = "magma", limits = c(0, 0.7)) +
  #scale_fill_gradient2(
    #low = "darkblue", mid = "white", high = "red", midpoint = 0.3, limits = c(0, 0.7) 
  #) +
  scale_fill_gradientn(
    colours = c("#000066", "#1919FF", "gold", "red"),
    values = scales::rescale(c(0, 0.2, 0.34, 0.7)),
    limits = c(0, 0.7)
  ) +
  scale_alpha_manual(
    values = c("ns" = 0.4, "significant" = 1)
  ) +
  # geom_text(aes(label = marker),
  #           color = "black",
  #           size = 7,           
  #           fontface = "bold",
  #           na.rm = TRUE) +
  # Dot for 0.3 < effsize <= 0.5
  geom_point(
    data = subset(all_wilcox_df, effsize > 0.3 & effsize <= 0.5),
    shape = 16, size = 2, color = "black",
    show.legend = FALSE
  ) +
  
  # Star for effsize > 0.5
  geom_point(
    data = subset(all_wilcox_df, effsize > 0.5),
    shape = 8, size = 2, color = "black",
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(
    title = ext,
    x = "CNV",
    y = "pathway",
    fill = "effect size",
    alpha = "stat. significance (p < 0.05)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    #panel.grid = element_blank()
    panel.grid.major = element_line(color = "grey90", linewidth = 0.05)
    )

ggsave(
  paste0(outdir, "/PDF_NEW_02_wilcoxtest_CNVxPW_v2.pdf"),
  p,
  width = 12.5, 
  height = 6, 
  dpi = 400)

# for further downstream analysis:
# identify cnvs of interest - cnvs that were associated with the activity of a certain pathway 
# (significant p value and moderate effect size)
# SINCE 0.3 is harsh cutoff, take cutoff value slightly below (effect sizes of 0.295 vs. 0.301 might be just as interesting ...)
# cnvs with moderate effect size as df
# 2025-11-06: change cutoff to 0.3
mod = all_wilcox_df %>%
  filter(effsize > 0.3) 

# get cnv id (to subset "cnv_list", keep infos to identify which barcodes had a given cnv vs. which did not)
vals = unique(mod$cnv)
sublist = cnv_list[vals]
wilcox_sublist = wilcox_list[vals]

# save "mod" df and sublist
#write.csv(
  #mod, 
  #paste0(outdir, "/NEW_03_overview_CNVxPW.csv"), 
  #row.names = FALSE)

#qsave(
  #sublist, 
  #paste0(outdir, "/NEW_04_bc_info_cnvs_assoc_w_pws.qs")
#)

# also save wilcoxon test results, just to be safe
#qsave(
  #wilcox_sublist, 
  #paste0(outdir, "/NEW_05_wilcox_test_CNVxPW.qs")
#)
