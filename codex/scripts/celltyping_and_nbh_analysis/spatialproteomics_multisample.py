#########################################################################
# CODEX thresholding, binarization, celltyping and neighborhood analysis
#########################################################################

# load packages
import spatialproteomics as sp
from skimage.io import imread
import matplotlib.pyplot as plt
from scipy.signal import medfilt2d
import pandas as pd
import time
import yaml
import xarray as xr
import numpy as np
import os

##############
# load images
##############

# imgs were kept in vals folder
dir = "/g/saka/Valeriy/projects/tFL/codex/workflow/intermediate/"

# set outdir for generated figs
figdir = "/g/saka/Tatjana/data/04_CODEX/figs/final/nbh_composition/"
# make separate dir fr saving the figs of celltypes and neighborhoods identified
figdir_ct ="/g/saka/Tatjana/data/04_CODEX/figs/final/ct/"
figdir_nbh = "/g/saka/Tatjana/data/04_CODEX/figs/final/nbh/"

# load list w names of zarr files
with open("/g/saka/Tatjana/analysis/codex/codex_roinames_LN0027_LN0438.txt", "r") as f:
    rois = f.read().splitlines()

# make path to imgs out of base dir and roi names 
zarr_pat = [os.path.join(dir, f"{name}.zarr") for name in rois]
print(zarr_pat)

# load all imgs
im = [xr.open_dataset(pat, engine = "zarr") for pat in zarr_pat]

zarr_pat
im

# subset to test script
# FOR RUNNING SCRIPT INTERACTIVELY:
# subset im to only contain 1-2 samples
# otherwise will take realy long

#im = im[2:4]
#rois = rois[2:4]



##############################
# add cell segmentation masks
##############################

## NOTE: done outside of this script, therefore outcomment code snippet! Same for thresholding,
## was applied outside of this script so is outcommented here. NOTE: polish scripts and add to folder!!

# load segmentation masks 
#sm_pat = [os.path.join(dir, f"segmentation/merged/{name}.npy") for name in rois]
#sm = [np.load(pat) for pat in sm_pat]

#sm_test = sm[0:3]

# add segmentation masks to images
# im_test = [img.pp.add_segmentation(seg) for img, seg in zip(im_test, sm_test)]

##############
# crop images
##############
# otherwise run out of memory when wanting to do things interactively.

#crop_im = [
    #img.pp[5000:8000,4000:7000]
    #for img in im
#]

###################
# img thresholding
###################

# for thresholding imgs: load thresholding yaml
# suitable thresholds were manually determined for each img using napari and added to yml
# with open("/g/saka/Tatjana/analysis/codex/yamlfiles/thresholds.yaml") as f:
#     thr =  yaml.safe_load(f)

# subset thresholds by roi of interest - since here, we only have 1 roi = 1 img
# th = thr['thresholds']['LN0438_MAAFHY1_R3']
# th

# issue that some channels in "th" are empty, = they were excluded because of bad antibody quality.
# fill up empty strings with "1" (="threshold out 100%", = marker will be excluded from analysis)
# for k, val in th.items():
#     if val == '':
#         th[k] = 1


# # now apply thresholds to images
# th_im = [
#     img.pp.threshold(
#         quantile=[thr[ch] for ch in img['channels'].data],
#         channels=img['channels'].data,
#         #key_added='_thresh_img',
#         shift=False
#     )
#     for img, thr in zip(im, th)
# ]

#################################
# quantifying protein expression
#################################

# use segmentation masks to quantify protein expression
# taking median intensity and applying arcsinh transform 
# see spatialproteomics vignette:
# https://sagar87.github.io/spatialproteomics/notebooks/CellTypePrediction.html

im = [
    img.pp.add_quantification(func="intensity_mean") 
    .pp.transform_expression_matrix(method="arcsinh")
    for img in im
]

########################################
# broad celltype prediction using argmax
########################################

# make preliminary celltyping dict: basic celltypes, no subtypes
# note: excluded two markers that were very dominant in img, probably a thresholding issue -> revisit
mk_dict = {
    'B': "PAX5",
    'B': "CD20",
    'T': "CD3",
    'Myeloid': "CD11b",
    'Dendritic': "CD11c",
    'Granulo': "CD15",
    'Macro': "CD68",
    'NK': "CD56",
    'FDC': "CD35",
    #'Endo': "CD31",
    'FRC': "CD90"
    #'B_plasma': "CD138"
}

# celltype prediction
im = [
    img.la.predict_cell_types_argmax( 
        mk_dict, key = "_intensity")
    # add colors that we had specified above
    #.la.set_label_colors(list(mk_dict.keys()), colors)
    for img in im
]

## addition 28/10/25: want to export images with single marker stainings of some
## markers that are super interesting, but they are very dim: 
# apply gamma correction to imgs to make them brighter

'''
# make plots for marker gene expression "PAX5", "CD3", "CD11c", "CD90"
ch = ["PAX5", "CD3", "CD11c", "CD90", "CD20"]

figdir2 = "/g/saka/Tatjana/data/04_CODEX/figs/single_marker_expressions_white_gamma0.05/"

for i, roi in enumerate(im):
    for c in ch:
        img = roi.pp[c]["_image"]
        #print(img)
        #img_norm = (img - img.min()) / (img.max() - img.min())
        gamma =0.05
        img_gamma = img ** gamma
        img_gamma = img_gamma[0, :, :]
        
        plt.figure(figsize = (7, 7))
        plt.imshow(img_gamma, cmap="gray")
        plt.axis("off")
        #_ = roi.pp[c].pl.colorize("#FDFEFF").pl.show()

        plt.title(f"sample {i}, channel {c}")
        savepat = os.path.join(figdir2, f"{i}_{c}.png")

        plt.savefig(savepat, dpi = 300, bbox_inches = "tight")
        plt.close()

'''


#########################################
# marker binarization for cell subtyping
#########################################

# above, we defined basic cell types (B cells, T cells, Macrophages, ...)
# now, go one level deeper and predict cell subtypes: 
# for this, necessary to binarize markers (i.e., is cell pos. or neg. for that given marker?)
# done with quantiles: e.g. quantile 0.7 = more than 70% of pixels of a given cell must be positive for 
# a given marker, only then cell counts as "positive" 
# tricky with cell surface markers: since most imgs of cells are cross-sections, often only the outer rim of a cell
# is pos for cell surface markers: in this case, threshold was set very low (e.g. 0.1, 0.15 = 10%/15% positive pixels 
# necessary for cell to be "positive")
# performed binarization manually in napari, load binarization quantiles as yaml file

# NOTE: for running interactively, used subset of images. 
# Also use subset of rois (for accessing correct binarization thresholds 
# rois_s = rois[0:2]

with open("/g/saka/Tatjana/analysis/codex/yamlfiles/binarization.yaml") as f:
    bina =  yaml.safe_load(f)
# subset bin to sample that we are processing, and copy it (since we are applying the same binarization
# quartiles to both our mock samples)

bn = [bina['thresholds'][roi] for roi in rois]
bn

# issue that some channels in "bn" are empty, = they were excluded, fill up empty strings with "1"
# (="threshold out 100%", = marker will be excluded from analysis)
bn = [{k: (1 if v == '' else v) for k, v in dct.items()} for dct in bn]

# now binarize functional markers, done in multistep process:
# 1.) again quantify protein expression (instead of mean expression use % positive pixels now)
# 2.) perform binarization with la.threshold_labels()
im = [img.pp.add_quantification(
    func = sp.percentage_positive, key_added = "_percentage_positive")
         for img, b in zip(im, bn)]

im = [
    img.la.threshold_labels(b, layer_key = "_percentage_positive")
    for img, b in zip(im, bn)
]

# inspect binarization results
# im_ct_bin[0].pp.get_layer_as_df().head()


#################
# cell subtyping 
#################

# combine basic cell type predictions with marker binarizations to predict cell subtypes
# decision tree to define cell type hierarchy, imported as yaml file
with open("/g/saka/Tatjana/analysis/codex/yamlfiles/celltyping_final.yml") as f:
    ct =  yaml.safe_load(f)

# make color dict for cell subtyping
type_cols = {
    'B': "#ad1d48",
    'T': "#060d70",
    'Granulo': "#543403",
    'Macro': "#7A33FF",
    'NK': "#E849A5",
    'FDC': "#0000FF", 
    #'Endo': "#2FA12F"
    'FRC': "#8af293", 
    'Myeloid': "#D98A00", 
    'Dendritic': "#FFB000",
    #'B_plasma': "#9E66FF",
    'T_h': "#11B6CE",        
    'T_h_naive': "#fb9a99",   
    'T_h_mem': "#a6cee3",      
    'T_reg_CD25+FoxP3+': "#b2df8a",  
    'T_reg_CD25+FoxP3-': "#cab2d6",  
    'T_reg_CD25-FoxP3+': "#fdbf6f",  
    'T_tox': "#1b9e77",        
    'T_tox_naive': "#9955BB", 
    'T_tox_mem': "#7570b3"     
}



# add cell subtypes based on celltype yaml above
im = [
    img.la.predict_cell_subtypes(ct)
    .la.set_label_colors(type_cols.keys(), type_cols.values())
    for img in im
]


# plotting the cell type predictions
plt.figure(figsize=(10, 10))
_ = th_im_ct_bin[0].pl.autocrop().pl.show(render_image=False, render_labels=True)


 for i, roi in enumerate(im):
     plt.figure(figsize = (10, 10))
     _ = roi.pl.render_labels()
     _ = roi.pl.show(render_image=False, render_labels=True)

     ax = plt.gca()
     leg = ax.get_legend()
     if leg is not None:
         #leg.set_title(leg.get_title().get_text(), prop={'size': 8})
         #for text in leg.get_texts():
             #text.set_fontsize(8) 
         leg.remove() # remove legend for 1 trial run, uncomment this normally 

    
     plt.title(f"Cell type predictions {i}")
     savepat = os.path.join(figdir, f"{i}_nolegend.png")

     plt.savefig(savepat, dpi = 300, bbox_inches = "tight")
     plt.close()
    
    

#plt.savefig("/g/saka/Tatjana/analysis/codex/scripts/test_fig.png", dpi=300, bbox_inches="tight")
#plt.close()


########################
# neighborhood analysis
########################

# set neighborhood colors
neighborhoods, colors = [f"Neighborhood {x}" for x in np.arange(4)], ["#e41a1c", "#4daf4a", "#377eb8", "#984ea3"]

#im_test = im[0].nh.compute_neighborhoods_knn

# make dict of imgs
sp_dict = {i + 1: val for i, val in enumerate(im)}

# make ImageContainer
img_cont = sp.ImageContainer(sp_dict)

# calculate neieghborhoods
img_cont = img_cont.compute_neighborhoods(neighborhood_method='radius', radius=100, k = 4)
#img_cont = sp.ImageContainer(img_cont)


# make figs of nbh for all rois
# for key, val in sp_dict.items():
#     plt.figure(figsize = (10, 10))
#     _ = val.pl.autocrop().pl.show(render_image = False, render_neighborhoods=True)
    #savepat = os.path.join(figdir_nbh, f"{key}.png")
    
    #plt.savefig(savepat, dpi = 300, bbox_inches = "tight")
    #plt.close()
    

#.nh.set_neighborhood_colors(neighborhoods, colors)


# for each image in im:

#obs2 = img_cont[2].pp.get_layer_as_df("_obs")
#nbh_df2 = pd.crosstab(obs2['_labels'], obs2['_neighborhoods'])

for img, roi_name in zip(img_cont, rois):

    obs = img_cont[img].pp.get_layer_as_df("_obs")
    nbh_df = pd.crosstab(obs['_labels'], obs['_neighborhoods'])

    outpat = os.path.join(figdir, f"{roi_name}.csv")
    
    nbh_df.to_csv(outpat, index=True)






