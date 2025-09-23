#!/bin/bash
#SBATCH --job-name=inferCNV_array
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tatjana.tomek@embl.de
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=10-00:00:00
#SBATCH --array=0-21
#SBATCH --error=.slurm/%A_%a_siCNV.err
#SBATCH --output=.slurm/%A_%a_siCNV.out

module load JAGS
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
source /home/tomek/.Rprofile

# get roi_id (extension) from all_rois_list.txt
EXT=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" //g/saka/Tatjana/data/roinames_visium_prelimgrowth2.txt)

DATA_DIR="/scratch/tomek/01_siCNV/infercnv_outs_spotbased/NEW_spatial_doms_prelim_growth2/jobarray_all_pat_per_roi_spotbased_BMPN_0.3"
CNV_regs="${DATA_DIR}/${EXT}/17_HMM_predHMMi6.hmm_mode-cells.pred_cnv_regions.dat"
MCMC_obj="${DATA_DIR}/${EXT}/BayesNetOutput.HMMi6.hmm_mode-cells/MCMC_inferCNV_obj.rds"
OUTDIR="/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2/"

echo "Running job for $EXT"
Rscript /g/saka/Tatjana/analysis/visium/01_CNV_calling/infercnv_visium_spot-based/SPOT_01_infercnv_output_reformating.R "$CNV_regs" "$MCMC_obj" "$OUTDIR" "$EXT"
