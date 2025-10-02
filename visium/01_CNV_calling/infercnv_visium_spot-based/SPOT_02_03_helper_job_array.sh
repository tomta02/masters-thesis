#!/bin/bash
#SBATCH --job-name=find_common_cnvs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tatjana.tomek@embl.de
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=01:00:00
#SBATCH --array=0-21
#SBATCH --error=.slurm/%A_%a_siCNV.err
#SBATCH --output=.slurm/%A_%a_siCNV.out

module load JAGS
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
source /home/tomek/.Rprofile

# get roi_id (extension) from all_rois_list.txt
EXT=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" //g/saka/Tatjana/data/roinames_visium_prelimgrowth2.txt)

DATA_DIR="/g/saka/Tatjana/data/01_CNV_analysis/01_infercnv_output/spotbased/all_pat_analyzed_spotbased_BMPN_0.3_NEW_prelim_growth2"
CNV_regs="${DATA_DIR}/${EXT}/SPOT_01_cnv_final_cnvprob_bigger_than_0.9.csv"

echo "Running job for $EXT"
Rscript /g/saka/Tatjana/analysis/visium/01_CNV_calling/infercnv_visium_spot-based/SPOT_02_find_common_cnvs_w_genomic_ranges.R "$DATA_DIR" "$CNV_regs" "$EXT"
