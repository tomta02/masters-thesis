#============================================================================#
# Decontamination of Visium spots using SpotClean      
#============================================================================#

# AIM: using a healthy lymph node dataset from 10X as healthy control dataset
# for comparison to our follicular lymphoma samples. 
# Subjecting it to the same prepro as patient datas (i.e. SpotClean + basic QC)
# here: SpotClean
# healthy LN dataset: https://www.10xgenomics.com/datasets/human-lymph-node-1-standard-1-1-0
# SpotClean: https://github.com/zijianni/SpotClean

# load libraries
library(SpotClean)
library(qs)
library(SpatialExperiment)


set.seed(99)


# test
# read in healthy lymph node spaceranger "outs" as seurat object
healthy_lymph_test <- Load10X_Spatial(
  data.dir='/g/saka/Tatjana/visium/V1_healthy_human_lymph_node/outs',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Visium"
)

rownames(healthy_lymph_test)[1:5]

# START
slide_obj <- read10xVisium(samples = "/g/saka/Tatjana/visium/V1_healthy_human_lymph_node/outs", 
                           data = "raw")
decont_obj <- spotclean(
  slide_obj,
  gene_cutoff = 0.1
  )

summary(metadata(decont_obj)$contamination_rate)


# convert slide object to Seurat object for downstream analysis
seu_obj <- convertToSeurat(
  decont_obj,
  image_dir = "/g/saka/Tatjana/visium/V1_healthy_human_lymph_node/outs/spatial",
  )

rownames(seu_obj)[1:5]

# save Seurat object
qsave(
  seu_obj, 
  "/g/saka/Tatjana/visium/V1_healthy_human_lymph_node/healthy_ln_spotcleaned_2.rds"
  )
