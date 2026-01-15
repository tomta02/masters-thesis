#################################
# CNVxPW_02_downstream_analysis
#################################

# calculating a null-distribution for the CNV x pathway associations that were
# found interesting. I.e.: we found that the Visium spots that were having a 
# certain CNV were having a significantly different pathway activation score 
# for certain pathways compared to the spots that did not have a certain CNV
# plot this out: i.e 
# 1.) compare PW activation levels of spots WITH and WITHOUT CNV
# 2.) calculate a null distribution: i.e. if our sample had e.g. 1000 spots, 
# and in 800 of them there is a certain CNV present and in 200 not, take 1000x
# 800 random spots vs. 200 random spots and calculate pathway activation level
# differences between them, and see whether pw activation values differ

args = commandArgs(trailingOnly = TRUE)

dir = args[1] # "/g/saka/Tatjana/data/03_CNVxPW_analysis/" #
ext = args[4] # "LN0438_MAAFHY3_R2" # 
cnvs_mod_effsize = args[2] # paste0(dir, ext, "/NEW_03_overview_CNVxPW.csv") #
pw_bc = args[3] # paste0(dir, ext, "/NEW_04_bc_info_cnvs_assoc_w_pws.qs") #



# load libraries
library(qs)
library(tidyverse)
library(ggplot2)
library(rstatix)

# 1. load output data of CNVxPW_02

# load CNVs that were associated with at least a moderate effect size difference
# i.e. moderately different pathway activation levels between spots that had a givne
# cnv and spots that did not
mod = read.csv(cnvs_mod_effsize)


# also load barcode-level information
bc_pw_info = qread(
  pw_bc,
  nthreads = 32
)


# add identifying cnv column (pattern: "chrX_start_stop")
for (i in seq_along(bc_pw_info)) {
  df = bc_pw_info[[i]]
  df$cnv = paste(df$chr, (round(df$start/1000)), (round(df$end/1000)), sep = "_")
  bc_pw_info[[i]] = df
}


# get pathways that are interesting
pws_tocheck = unique(mod$variable)

pw_results = map(bc_pw_info, ~
                    map_dfr(pws_tocheck, function(v) {
                      tmp = .x %>%
                        group_by(cnv_test) %>%
                        summarise(mean = mean(.data[[v]], na.rm = TRUE),
                                  sd = sd(.data[[v]], na.rm = TRUE), 
                                  n = sum(!is.na(.data[[v]])),
                                  .groups = "drop")
                    
                      
                      # for each df, need to extract chromosom "chr", "start', "stop" values
                      # to later identify the CNVs when plotting
                      # are the same across all rows, so just do unique() 
                      chr_v = unique(.x$chr)
                      start_v = unique(.x$start)
                      stop_v = unique(.x$end)
                      
                      pval = wilcox.test(
                        na.omit(.x[[v]][.x$cnv_test == 0]),
                        na.omit(.x[[v]][.x$cnv_test == 1]),
                      )$p.value
  
                      tibble(
                        variable = v,
                        cnv = as.character(paste(chr_v, (round(start_v/1000)), (round(stop_v/1000)), sep ="_")),
                        group0_mean = tmp$mean[tmp$cnv_test == 0],
                        group0_sd = tmp$sd[tmp$cnv_test == 0],
                        group0_n = tmp$n[tmp$cnv_test == 0],
                        group1_mean = tmp$mean[tmp$cnv_test == 1],
                        group1_sd = tmp$sd[tmp$cnv_test == 1],
                        group1_n = tmp$n[tmp$cnv_test == 1],
                        mean_diff = group1_mean - group0_mean,
                        pval = pval
                      )
                    })
)


# test = bc_pw_info_1[bc_pw_info_1$cnv_test == 1,]
# mean(test$NFkB, na.rm = TRUE)


# Helper function to generate plots for a given variable
plotPW = function(dflist, results_list, var) {
  imap(dflist, function(df, name) {
    # extract summary stats
    stats = results_list[[name]] %>% filter(variable == var)

    # keep only CNV 0/1 and non-missing values for this variable
    df = df %>% filter(cnv_test %in% c(0, 1)) %>%
      filter(!is.na(.data[[var]]))

    cnv_title = df$cnv


    ggplot(df, aes(x = factor(cnv_test), y = .data[[var]], fill = factor(cnv_test))) +
      geom_violin(trim = FALSE, alpha = 0.5) +
      geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
      # overlay means and SDs
      geom_point(data = stats, aes(x = 1, y = group0_mean), color = "red", size = 3, inherit.aes = FALSE) +
      geom_point(data = stats, aes(x = 2, y = group1_mean), color = "red", size = 3, inherit.aes = FALSE) +
      geom_errorbar(data = stats,
                    aes(x = 1, ymin = group0_mean - group0_sd, ymax = group0_mean + group0_sd),
                    color = "red", width = 0.1, inherit.aes = FALSE) +
      geom_errorbar(data = stats,
                    aes(x = 2, ymin = group1_mean - group1_sd, ymax = group1_mean + group1_sd),
                    color = "red", width = 0.1, inherit.aes = FALSE) +
      scale_fill_manual(values = c("0" = "darkgray", "1" = "orange")) +
      labs(title = paste(var, "PW activity by CNV status:", cnv_title),
           x = "CNV absence: 0\nCNV_presence: 1",
           y = paste(var, "pathway activity value"),
           fill = "cnv_test") +
      theme_minimal(base_size = 14)
  })
}


# generate plots, store in list
p_l = list()

for (pw in pws_tocheck) {
  p_l[[pw]] = plotPW(bc_pw_info, pw_results, pw)
}

# now save plots
# mayve think how to rewrite tidyverse style
# loop over all pathways, get list of ggplots
for (pwname in names(p_l)) {
  pw_p = p_l[[pwname]]

  # save each plot individuallt
  for (i in seq_along(pw_p)) {

    fname = paste0(dir, ext, "/figs/violinplot_",
                       pwname, "_", i, ".png")

    ggsave(
      filename = fname,
      plot = pw_p[[i]],
      width = 9,
      height = 4.5,
      dpi = 400
    )
  }
}



####################################################################################
# next: calculate null distribution
# see how likely our observed effect size is to get when shuffling group assignments
# (i.e. cnv_test = 1 vs. cnv_test = 0; i.e. whether a given barcode has a CNV 
# or not) while keeping group sizs constant. = permutation test
####################################################################################

### IMPORTANT: RERUN WITH JAK STAT LATER
pws_tocheck2 = pws_tocheck[pws_tocheck != "JAK-STAT"]

set.seed(555)

# nr of permutations, start w 1000
nperm = 3000
permresults = list()

for (i in seq_along(bc_pw_info)) {
  df = bc_pw_info[[i]] %>% filter(cnv_test %in% c(0,1))
  
  chr_v = unique(df$chr)
  start_v = unique(df$start)
  end_v = unique(df$end)
  cnv_v = unique(df$cnv)
  
  
  permresults[[i]] =lapply(pws_tocheck2, function(var) {
    
    # remove rows where the variable currently being tested is NA
    df2 = df %>% filter(!is.na(.data[[var]]))
    
    # calculate group sizes after removing NAs
    n0 = sum(df2$cnv_test == 0)
    n1 = sum(df2$cnv_test == 1)
    
    # get effect size of the actual group asignments i.e. cnv_test = 0 vs. cnv_test = 1
    obs_eff = wilcox_effsize(df2,
                             as.formula(paste(var, "~ cnv_test")), 
                             ci = FALSE)$effsize
    
    # get usable values
    allvals = df2[[var]]
    
    # get permutation null distribution
    null_eff = replicate(nperm, {
      shuffled = sample(allvals)
      group = c(rep(0, n0), rep(1, n1))
      wilcox_effsize(
        data.frame(value = shuffled, group = group),
        value ~ group,
        ci = FALSE
      )$effsize
    })
    
    # last, calc permutation pval
    #perm_p = (1 + sum(abs(null_eff) >= abs(obs_eff))) / (1 + nperm)
    perm_p = (1 + sum(abs(null_eff) >= abs(obs_eff))) / (1 + nperm)
    
    # print(var)
    # print(chr_v)
    # print(start_v)
    # print(end_v)
    # print(cnv_v)
    # print(perm_p)
    # print(mean(abs(null_eff)))
    
    data.frame(
      variable = var,
      chr = chr_v,
      start = start_v,
      end = end_v,
      cnv = cnv_v,
      observed_eff = obs_eff,
      perm_p = perm_p,
      mean_null_eff = mean(abs(null_eff))
      )

  }) %>% bind_rows()

}

# combine res
permdf = bind_rows(permresults)

# save permdf
write_csv(permdf, paste0(dir, ext, "/NEW_permutation_test_results_20251128.csv"))

# also save pw_res_df
pw_results_df = do.call("rbind", pw_results)
write_csv(pw_results_df, paste0(dir, ext, "/NEW_pw_stats_20251128.csv"))



