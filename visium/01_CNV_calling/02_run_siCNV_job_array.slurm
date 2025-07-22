#!/bin/bash
#SBATCH --job-name=inferCNV_array
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tatjana.tomek@embl.de
#SBATCH --partition=bigmem
#SBATCH --constraint=turin
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=10-00:00:00
#SBATCH --array=0-21
#SBATCH --error=.slurmR/%A_%a_siCNV.err
#SBATCH --output=.slurmR/%A_%a_siCNV.out

module load JAGS
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
source /home/tomek/.Rprofile

# get roi_id (extension) from all_rois_list.txt
EXT=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" /g/saka/Tatjana/visium/analysis_1/01_CNV_analysis/patient_data_formated_for_running_infercnv/all_pat_countmat_and_annot_split_by_roi/all_rois_list.txt)

DATA_DIR="/g/saka/Tatjana/visium/analysis_1/01_CNV_analysis/patient_data_formated_for_running_infercnv/all_pat_countmat_and_annot_split_by_roi"
GENEORDER_F="/g/saka/Tatjana/visium/analysis_1/01_CNV_analysis/hg38_gencode_v27.txt"
RAW_COUNTS="${DATA_DIR}/${EXT}_countmtx.tsv"
ANNOTS="${DATA_DIR}/${EXT}_annot.tsv"
PATIENT_ID=$(echo "$EXT" | cut -d'_' -f1)
OUTDIR="/scratch/tomek/01_siCNV/RESULTS_spotbased_infercnv_per_roi/jobarray_all_pat_roi_per_spot_BayesMaxPNormal_0.5/"
SCRIPT_PATH="/g/saka/Tatjana/visium/analysis_1/scripts/01_CNV_calling/02_run_siCNV.R"

echo "Running job for $EXT"
Rscript "$SCRIPT_PATH" "$RAW_COUNTS" "$GENEORDER_F" "$ANNOTS" "$OUTDIR" "$EXT"
