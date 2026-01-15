#!/bin/bash
#SBATCH --job-name=siCNV
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tatjana.tomek@embl.de
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=32G
#SBATCH --time=10-00:00:00
#SBATCH --error=.slurm/%j_siCNV.err
#SBATCH --output=.slurm/%j_siCNV.out
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
source /home/tomek/.Rprofile
Rscript -e "
  rds_file = readRDS('/scratch/tomek/01_siCNV/RESULTS_spotbased_infercnv_per_roi/jobarray_all_pat_roi_per_spot_BayesMaxPNormal_0.5/LN0025_1DA209_R1/BayesNetOutput.HMMi6.hmm_mode-cells/MCMC_inferCNV_obj.rds')
  saveRDS(rds_file@cell_probabilities, file='/scratch/tomek/01_siCNV/RESULTS_spotbased_infercnv_per_roi/jobarray_all_pat_roi_per_spot_BayesMaxPNormal_0.5/LN0025_1DA209_R1/cellprobs_mcmc.rds')
  saveRDS((colnames(rds_file@count.data)), file='/scratch/tomek/01_siCNV/RESULTS_spotbased_infercnv_per_roi/jobarray_all_pat_roi_per_spot_BayesMaxPNormal_0.5/LN0025_1DA209_R1/bc_names.rds')
  saveRDS(rds_file@cnv_regions, file='/scratch/tomek/01_siCNV/RESULTS_spotbased_infercnv_per_roi/jobarray_all_pat_roi_per_spot_BayesMaxPNormal_0.5/LN0025_1DA209_R1/cnv_regions.rds')
  saveRDS(rds_file@cell_gene, file='/scratch/tomek/01_siCNV/RESULTS_spotbased_infercnv_per_roi/jobarray_all_pat_roi_per_spot_BayesMaxPNormal_0.5/LN0025_1DA209_R1/cell_gene.rds')
"
