#!/bin/bash
#SBATCH --job-name=siCNV_roibased
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tatjana.tomek@embl.de
#SBATCH --partition=bigmem
#SBATCH --constraint=turin
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=832G
#SBATCH --time=10-00:00:00
#SBATCH --error=.slurm/%j_siCNV.err
#SBATCH --output=.slurm/%j_siCNV.out
RAW_COUNTS="/g/saka/Tatjana/data/01_CNV_analysis/patient_data_formated_for_running_infercnv/all_pat_countmat_and_annot_new_spatial_domains/LN0438_countmtx.tsv"
GENEORDER_F="/g/saka/Tatjana/data/01_CNV_analysis/hg38_gencode_v27.txt"
ANNOTS="/g/saka/Tatjana/data/01_CNV_analysis/patient_data_formated_for_running_infercnv/all_pat_countmat_and_annot_new_spatial_domains/LN0438_annot.tsv"
OUTDIR="/scratch/tomek/01_siCNV/RESULTS_spatial_domain_based_new_prelim_growth2_BMPN_0.3/"
EXT="LN0438"

module load JAGS
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
source /home/tomek/.Rprofile
Rscript /g/saka/Tatjana/analysis/visium/01_CNV_calling/02_run_siCNV.R "$RAW_COUNTS" "$GENEORDER_F" "$ANNOTS" "$OUTDIR" "$EXT"