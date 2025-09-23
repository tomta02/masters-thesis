#!/bin/bash
#SBATCH --job-name=spotbased_infercnv_output_reformat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tatjana.tomek@embl.de
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --error=.slurm/%j_siCNV.err
#SBATCH --output=.slurm/%j_siCNV.out
CNV_REGS_PRELIM="/scratch/tomek/01_siCNV/infercnv_outs_spotbased/NEW_spatial_doms_prelim_growth2/jobarray_all_pat_per_roi_spotbased_BMPN_0.3/LN0438_MAAFHY1_R1/17_HMM_predHMMi6.hmm_mode-cells.pred_cnv_regions.dat"
MCMC_obj="/scratch/tomek/01_siCNV/infercnv_outs_spotbased/NEW_spatial_doms_prelim_growth2/jobarray_all_pat_per_roi_spotbased_BMPN_0.3/LN0438_MAAFHY1_R1/BayesNetOutput.HMMi6.hmm_mode-cells/MCMC_inferCNV_obj.rds"
OUTDIR="/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2/"
EXT="LN0438_MAAFHY1_R1/"

module load JAGS
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
source /home/tomek/.Rprofile
Rscript /g/saka/Tatjana/analysis/visium/01_CNV_calling/infercnv_visium_spot-based/SPOT_01_infercnv_output_reformating.R "$CNV_REGS_PRELIM" "$MCMC_obj" "$OUTDIR" "$EXT"
