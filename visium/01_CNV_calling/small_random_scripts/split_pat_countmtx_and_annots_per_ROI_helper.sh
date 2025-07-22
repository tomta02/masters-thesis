#!/bin/bash
VISLIST="/path/to/list/of/visium/seurat/objects/02_10_visium_flt.qs"
ANNOT="/path/to/patient_data_formated_for_running_infercnv/all_pat_countmat_and_annot/pat1_annot.tsv"
COUNTMTX="/path/to/patient_data_formated_for_running_infercnv/all_pat_countmat_and_annot/pat1_countmtx.tsv"
OUTDIR="/path/to/out/dir/all_pat_countmat_and_annot_split_by_roi/"

module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
source /home/tomek/.Rprofile
Rscript /path/to/rscripts/split_pat_countmtx_and_annots_per_ROI.R "$VISLIST" "$ANNOT" "$COUNTMTX" "$OUTDIR"
