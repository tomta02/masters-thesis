#!/bin/bash
#SBATCH --job-name=scatterplot_maligb_cnv
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tatjana.tomek@embl.de
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=01:00:00
#SBATCH --array=0-21
#SBATCH --error=.slurm/%A_%a_siCNV.err
#SBATCH --output=.slurm/%A_%a_siCNV.out

module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
source /home/tomek/.Rprofile

# get roi_id (extension) from all_rois_list.txt
EXT=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" //g/saka/Tatjana/data/roinames_visium_prelimgrowth2.txt)

DATA_DIR="/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2/roi_based_cnvcalling"
CNV_clean="${DATA_DIR}/${EXT}/NEW_SPOT_03_cnvs_cleaned_75th_quantile_binarized.csv"
CNV_clean_bin="${DATA_DIR}/${EXT}/NEW_SPOT_06_cleaned_cnvs_per_bc_from_NEW_SPOT_03_cnvs_cleaned_binarized_for_mking_charmat.csv"

echo "Running job for $EXT"
Rscript /g/saka/Tatjana/analysis/visium/01_CNV_calling/small_random_scripts/scatterplots_malig_B_cells_vs_CNVs_per_spot.R "$DATA_DIR" "$CNV_clean" "$CNV_clean_bin" "$EXT"
