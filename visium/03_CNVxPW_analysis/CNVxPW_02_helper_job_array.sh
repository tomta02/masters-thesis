#!/bin/bash
#SBATCH --job-name=CNV
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tatjana.tomek@embl.de
#SBATCH --ntasks=1
#SBATCH --mem=400G
#SBATCH --time=10:00:00
#SBATCH --array=0-21
#SBATCH --error=.slurm/%A_%a_siCNV.err
#SBATCH --output=.slurm/%A_%a_siCNV.out

module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
source /home/tomek/.Rprofile

# get roi_id (extension) from all_rois_list.txt
EXT=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" //g/saka/Tatjana/data/roinames_visium_prelimgrowth2.txt)

DIR="/g/saka/Tatjana/data/03_CNVxPW_analysis/"
MOD="${DIR}/${EXT}/NEW_03_overview_CNVxPW.csv"
PW_BC="${DIR}/${EXT}/NEW_04_bc_info_cnvs_assoc_w_pws.qs"

echo "running job for $EXT"
Rscript /g/saka/Tatjana/analysis/visium/03_CNVxPW_analysis/CNVxPW_02_downstream_analysis_standardized.R "$DIR" "$MOD" "$PW_BC" "$EXT"
