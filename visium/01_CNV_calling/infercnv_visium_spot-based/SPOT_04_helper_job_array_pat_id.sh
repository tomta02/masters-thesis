#!/bin/bash
#SBATCH --job-name=make_charmat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tatjana.tomek@embl.de
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=480G
#SBATCH --time=95:00:00
#SBATCH --array=0-5
#SBATCH --error=.slurm/%A_%a_siCNV.err
#SBATCH --output=.slurm/%A_%a_siCNV.out

module load JAGS
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
source /home/tomek/.Rprofile

# get roi_id (extension) from all_rois_list.txt
EXT=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" //g/saka/Tatjana/data/pat_ids.txt)

DATA_DIR="/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2/whole_patient_cnvcalling/20251004/03_CNVs_GRANGES_CLEANED"
CNV_regs=$(ls "${DATA_DIR}"/"${EXT}_SPOT_03_cnvs_cleaned_50th_quantile_binarized".csv | head -n 1)
OUT="/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2/whole_patient_cnvcalling/20251004/04_CHAR_MAT_FOR_PLOTS"

echo "Running job for $EXT"
Rscript /g/saka/Tatjana/analysis/visium/01_CNV_calling/infercnv_visium_spot-based/SPOT_04_data_wrangling_for_making_heatmaps.R "$DATA_DIR" "$CNV_regs" "$EXT" "$OUT"
