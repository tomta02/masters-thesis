### Master's thesis project: Identification of tumor-enriched pathways and neighborhoods in transformed follicular lymphoma

This repository contains work from my masters thesis project conducted in the lab of Sinem Saka at EMBL Heidelberg (December 2024 - November 2025)

In the course of the thesis project, I analyzed a patient cohort of six patients diagnosed with follicular lymphoma transformation. For each of the patients, we had access to the following data modalities:
* spatial transcriptomics data (10x Visium)
* multiplexed immunofluorescence images (Akoya Spatial Proteomics, i.e. "CODEX")
* joint single-cell RNA-seq and ATAC-seq data (10x multiome)

The main objectives of the project were to:

1. **identify malignant regions in tissue samples**
I inferred the presence of copy number variations (CNVs) from the 10x Visium data and compared them to CNVs identified in the 10x multiome data. Malignant tissue regions annotated by a pathologist based on morphology served as control. 

2. **investigate cancer-specific biological pathway activity**
Using the prior knowledge network *PROGENy* and statistical models from the package *decoupleR*, I investigated pathway activity in tumor vs. non-tumor regions. I further was able to identify relationships between difference in pathway activity levels and the presence of copy number variations.

3. **characterize cell types present in the cancer regions and tumor microenvironment** 
Using the CODEX spatial proteomics data, I identified cellular neighborhoods of different cell type compositions across patients. I also compared celltypes identified in the spatial proteomics data to the cell types present in the Visium data (which were inferred using the 10x multiome dataset and the package *cell2location*) in order to assess agreement of the two modalities.

**Note:** this repository is work in progress, the project is ongoing and I remain in the Saka lab on a visitor contract. The scripts used for the CODEX preprocessing and analysis are being refactored and will be uploaded shortly, and this repository will be extended with example data obtained during the thesis work. 
