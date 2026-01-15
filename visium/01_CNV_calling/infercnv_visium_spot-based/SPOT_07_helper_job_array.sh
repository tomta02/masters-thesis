#!/bin/bash
#SBATCH --job-name=infercnv_numbat_compare
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
NB="/g/saka/Kristy/projects/composite/analysis/rdata/03_03_numbat_list.qs"
DATA_DIR="/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2/roi_based_cnvcalling"
INFCNV="${DATA_DIR}/${EXT}/NEW_SPOT_06_REFERENCE_CNVs_cleaned_from_NEW_SPOT_03_cnvs_cleaned_binarized.csv"

echo "running job for $EXT"
Rscript /g/saka/Tatjana/analysis/visium/01_CNV_calling/infercnv_visium_spot-based/SPOT_07_comparison_of_spotbased_CNV_calls_to_numbat_results.R "$DATA_DIR" "$NB" "$INFCNV" "$EXT"
