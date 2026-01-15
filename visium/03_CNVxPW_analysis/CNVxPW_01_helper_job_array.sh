#!/bin/bash
#SBATCH --job-name=cnv_pw_correl
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tatjana.tomek@embl.de
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=10:00:00
#SBATCH --array=0-21
#SBATCH --error=.slurm/%A_%a_siCNV.err
#SBATCH --output=.slurm/%A_%a_siCNV.out

module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
source /home/tomek/.Rprofile

# get roi_id (extension) from all_rois_list.txt
EXT=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" //g/saka/Tatjana/data/roinames_visium_prelimgrowth2.txt)

PWACT="/g/saka/Tatjana/data/02_pathway_analysis/progeny_decoupler/sct_normalized/pr_custom/pw_act_vals_and_pvals/01_pw_act_vals_all_roi_flt_pval_smaller_0.05.qs"
CNVS="/g/saka/Tatjana/data/03_CNVxPW_analysis/${EXT}/${EXT}_01_list_of_all_cnvs_found_NEW_shouldbesameasold.qs"
OUT="/g/saka/Tatjana/data/03_CNVxPW_analysis/${EXT}"

echo "running job for $EXT"
Rscript /g/saka/Tatjana/analysis/visium/03_CNVxPW_analysis/CNVxPW_01_correlating_pw_activity_w_cnvs.R "$PWACT" "$CNVS" "$EXT" "$OUT"
