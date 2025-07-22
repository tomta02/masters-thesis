# converting ENSEMBL IDs to gene names
# after running spotclean, our gene symbols are expressed as ENSEMBL IDs.
# for whatver reason
# use biomaRt to convert them back.

library(Seurat)
library(biomaRt)
library(dplyr)
library(qs)

# load seurat object
healthy_ln = qread("/g/saka/Tatjana/visium/V1_healthy_human_lymph_node/healthy_ln_spotcleaned2_QC.qs")
# remove SCT assay, since it is causing issues when subsetting seurat (see (1) a few lines below)
DefaultAssay(object = healthy_ln) <- "Spatial"
healthy_ln[["SCT"]] = NULL

# connect to ensembl database with biomart
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# extract the ensembl ids from seurat obj
ensembl_ids = rownames(healthy_ln)

# retrieve the mapping
gene_mapping = getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = ensembl_ids,
  mart = ensembl
)

# remove empty gene symbols, remove duplicaets
gene_mapping = gene_mapping %>% filter(hgnc_symbol != "")
gene_mapping = gene_mapping[!duplicated(gene_mapping$ensembl_gene_id), ]


# (1) subset Seurat object to only include matched ENSEMBL ids
healthy_ln_filtered = healthy_ln[rownames(healthy_ln) %in% gene_mapping$ensembl_gene_id, ]
# this is the error I get when not removing SCT assay:Error in .subscript.2ary(x, i, j, drop = TRUE) : subscript out of bounds

# Match and replace ENSEMBL IDs with gene names
new_rownames = gene_mapping$hgnc_symbol[match(rownames(healthy_ln_filtered), gene_mapping$ensembl_gene_id)]
rownames(healthy_ln_filtered) <- new_rownames

#> dim(healthy_ln@assays$Spatial)
#[1] 11965  4027
#> dim(healthy_ln_filtered@assays$Spatial)
#[1] 11722  4027

# now, our seurat object again has original HGNC gene names. Save it
qsave(healthy_ln_filtered, "/g/saka/Tatjana/visium/V1_healthy_human_lymph_node/healthy_ln_sc2_QC_HGNCgenenames.qs")

