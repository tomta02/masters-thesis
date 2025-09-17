# running inferCNV

args = commandArgs(trailingOnly = TRUE)

raw_counts_mtx = args[1]
geneorder_f = args[2]
annot_file = args[3]
outdir = args[4]
extension = args[5]

set.seed(99)

# load libraries
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(infercnv)
library(phylogram)
library(ape)
library(hdf5r)
library(SpatialInferCNV)
library(Matrix)

options("Seurat.object.assay.version" = "v3")
options(scipen = 100)

infCNV_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = raw_counts_mtx,
                                             gene_order_file = geneorder_f,
                                             annotations_file = annot_file,
                                             delim = "\t",
                                             ref_group_name = "GC")


infCNV_run <- infercnv::run(infCNV_obj,
                            cutoff = 0.1,
                            out_dir = paste0(outdir, extension),
                            BayesMaxPNormal = 0.3,
                            cluster_by_groups = FALSE, #needs to be set to T when running whole sample-based analysis instead of spot-based - otherwise will only report CNVs that are shared among all spatial domains
                            num_threads = 8,
                            analysis_mode = "cells",
                            #tumor_subcluster = TRUE,
                            HMM_report_by = "cell",
                            #tumor_subcluster_partition_method = "leiden",
                            HMM = TRUE,
                            denoise = TRUE)

