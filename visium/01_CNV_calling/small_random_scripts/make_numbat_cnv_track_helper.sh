#!/bin/bash
#SBATCH --job-name=plot_nb_track
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tatjana.tomek@embl.de
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=1:00:00
#SBATCH --array=0-5
#SBATCH --error=.slurm/%A_%a_siCNV.err
#SBATCH --output=.slurm/%A_%a_siCNV.out

module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
source /home/tomek/.Rprofile

# get patient id
PAT=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" //g/saka/Tatjana/data/pat_ids.txt)

NB="/g/saka/Kristy/projects/composite/analysis/rdata/03_03_numbat_list.qs"
GT="/g/saka/Kristy/projects/composite/analysis/rdata/00_numbat_gt_forTatjana.qs"

Rscript /g/saka/Tatjana/analysis/visium/01_CNV_calling/small_random_scripts/make_numbat_cnv_track.R "$NB" "$GT" "$PAT"
