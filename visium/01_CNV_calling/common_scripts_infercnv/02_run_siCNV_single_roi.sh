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
#SBATCH --error=.slurm/%A_%a_siCNV.err
#SBATCH --output=.slurm/%A_%a_siCNV.out

module load JAGS
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
source /home/tomek/.Rprofile

# get roi_id (extension) from all_rois_list.txt
EXT="LN0040HD_12BV2R_R3"

DATA_DIR="/g/saka/Tatjana/data/01_CNV_analysis/00_patient_data_infercnv_format/all_pat_countmat_and_annot_new_spatial_domains_prelim_growth2/split_by_roi"
GENEORDER_F="/g/saka/Tatjana/data/01_CNV_analysis/hg38_gencode_v27.txt"
RAW_COUNTS="${DATA_DIR}/${EXT}_countmtx.tsv"
ANNOTS="${DATA_DIR}/${EXT}_annot.tsv"
PATIENT_ID=$(echo "$EXT" | cut -d'_' -f1)
OUTDIR="/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/20251114_rerun_LN0040HD_12BV2R_R3/"
SCRIPT_PATH="/g/saka/Tatjana/analysis/visium/01_CNV_calling/common_scripts_infercnv/02_run_siCNV.R"

echo "Running job for $EXT"
Rscript "$SCRIPT_PATH" "$RAW_COUNTS" "$GENEORDER_F" "$ANNOTS" "$OUTDIR" "$EXT"
