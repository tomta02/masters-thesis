#!/bin/bash
#SBATCH --job-name=codex
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tatjana.tomek@embl.de
#SBATCH --partition=bigmem
#SBATCH --constraint=turin
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=800G
#SBATCH --time=93:00:00
#SBATCH --error=.slurm/%A_%a_siCNV.err
#SBATCH --output=.slurm/%A_%a_siCNV.out

source /g/korbel/tatjana_tomek/miniconda3/bin/activate

conda activate spatialproteomics
python /g/saka/Tatjana/analysis/codex/scripts/spatialproteomics_multisample.py



