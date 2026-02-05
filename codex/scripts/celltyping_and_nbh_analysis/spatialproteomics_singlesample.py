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
# load image
##############

img = xr.open_dataset("/g/saka/Tatjana/analysis/codex/zarr_file_example/LN0438_MAAFHY1_R3.zarr", engine="zarr")
img

##############################
# add cell segmentation masks
##############################

# load segmentation mask val had created
seg = np.load("/g/saka/Valeriy/projects/tFL/codex/workflow/segmentation/merged/LN0438_MAAFHY1_R3.npy")
seg.shape

# add segmentation masks to images
img = img.pp.add_segmentation(seg)


###################
# img thresholding
###################

# for thresholding imgs: load thresholding yaml
# suitable thresholds were manually determined for each img using napari and added to yml
with open("/g/saka/Tatjana/analysis/codex/yamlfiles/thresholds.yaml") as f:
    thr =  yaml.safe_load(f)

# subset thresholds by roi of interest - since here, we only have 1 roi = 1 img
th = thr['thresholds']['LN0438_MAAFHY1_R3']
th

# issue that some channels in "th" are empty, = they were excluded, fill up empty strings with "1"
# (="threshold out 100%", = marker will be excluded from analysis)
for k, val in th.items():
    if val == '':
        th[k] = 1
        
# threshold the img using quantiles defined in yaml
channels = img['channels'].data
thresholds = [th[ch] for ch in channels]

th_img = img.pp.threshold(
    quantile=thresholds,
    channels=channels,
    #key_added='thresh_img',
    shift=False)

#################################
# quantifying protein expression
#################################

# use segmentation masks to quantify protein expression
# taking median intensity and applying arcsinh transform 
# see spatialproteomics vignette:
# https://sagar87.github.io/spatialproteomics/notebooks/CellTypePrediction.html

th_img = th_img.pp.add_quantification(func="intensity_mean").pp.transform_expression_matrix(method="arcsinh")
th_img

# layer_key parameter not avail. Modify "_image" layer all the time - ask Val
#la = th_im.la

#fun = la.predict_cell_subtypes
#inspect.signature(fun)

#print(fun.__doc__)


########################################
# broad celltype prediction using argmax
########################################

# define channels of interest for now
# as well as colors
ch = ["PAX5", #B
      "CD3", #T
      "CD15", #Granulocyte
      "CD68", #Macrophage
      "CD56", #NK
      "CD35", #FDC
      "CD31", #Endothelial
      "CD90", #Fibrobalstic reticular cells
      "CD138"] #plasma cells

colors = [
    "#e6194B",
    "#3cb44b",
    "#ffe119",
    "#4363d8",
    "#ffd8b1",
    "#f58231",
    "#911eb4",
    "#fffac8",
    "#469990"]

# make preliminary celltyping dict: basic celltypes, no subtypes
mk_dict = {
    'B': "PAX5",
    'T': "CD3",
    'Granulo': "CD15",
    'Macro': "CD68",
    'NK': "CD56",
    'FDC': "CD35",
    'Endo': "CD31",
    'FRC': "CD90",
    'B_plasma': "CD138"
}

# celltype prediction, and adding colors that were specified above
th_img_ct = th_img.la.predict_cell_types_argmax(
    mk_dict, key = "_intensity").la.set_label_colors(list(mk_dict.keys()), colors)
    

# now plot celltype predictions next to markers
fig, ax = plt.subplots(1, 2, figsize=(10, 5))

_ = th_img_ct.pp[ch].pl.colorize(colors).pl.show(ax=ax[0])
_ = th_img_ct.pl.show(render_image=False, render_labels=True, ax=ax[1])


#########################################
# marker binarization for cell subtyping
#########################################

# above, we defined basic cell types (B cells, T cells, Macrophages, ...)
# now, go one level deeper and predict cell subtypes: 
# for this, we binarize markers (i.e., is cell pos. or neg. for that given marker?) NOTE: this is necessary bc the antibody staining in our imgs is crappy,
# in the original CODEX paper from Nolan they don't binarize

# done with quantiles: e.g. quantile 0.7 = more than 70% of pixels of a given cell must be positive for 
# a given marker, only then cell counts as "positive" 
# tricky with cell surface markers: since most imgs of cells are cross-sections, often only the outer rim of a cell
# is pos for cell surface markers: in this case, threshold was set very low (e.g. 0.1, 0.15 = 10%/15% positive pixels necessary for cell to be "positive")
# performed binarization manually in napari, load binarization quantiles as yaml file

with open("/g/saka/Tatjana/analysis/codex/yamlfiles/binarization.yaml") as f:
    bina =  yaml.safe_load(f)

bn = bina['thresholds']['LN0438_MAAFHY1_R3']
bn

# issue that some channels in "bn" are empty, = they were excluded, fill up empty strings with "1"
# (="threshold out 100%", = marker will be excluded from analysis)
for k, val in bn.items():
    if val == '':
        bn[k] = 1


# now binarize functional markers, done in multistep process:
# 1.) again quantify protein expression (instead of mean expression use % positive pixels now)
# 2.) perform binarization with la.threshold_labels()
th_img_ct = th_img_ct.pp.add_quantification(
    func = sp.percentage_positive, key_added = "_percentage_positive")

th_img_ct_bin = th_img_ct.la.threshold_labels(bn, layer_key = "_percentage_positive")

# inspect binarization results
th_img_ct_bin.pp.get_layer_as_df().head()


#################
# cell subtyping 
#################

# combine basic cell type predictions with marker binarizations to predict cell subtypes
# decision tree to define cell type hierarchy, imported as yaml file
with open("/g/saka/Tatjana/analysis/codex/yamlfiles/celltyping.yml") as f:
    ct =  yaml.safe_load(f)

# add cell subtypes based on celltype yaml above
th_img_ct_bin = th_img_ct_bin.la.predict_cell_subtypes(ct)



# plotting the cell type predictions
plt.figure(figsize = (10, 10))
_ = th_img_ct_bin.pl.autocrop().pl.show(render_image = False, render_labels = True)

plt.savefig("/g/saka/Tatjana/analysis/codex/scripts/test_fig.png", dpi = 300, bbox_inches = "tight")
plt.close()


########################
# neighborhood analysis
########################

# put our image(s) inside image container
# not sure I understand the purpose of this
# seems to be used for retrieving neighborhood composition across
# multiple samples, necessary for calling e.g. "get_neighborhood_composition".
# necessary to create dict of objects first
sp_dict = {"LN0438_MAAFHY1_R3": th_img_ct_bin}
img_cont = sp.ImageContainer(sp_dict)

img_cont = img_cont.compute_neighborhoods(k=5)


plt.figure(figsize = (10, 10))
_ = sp_dict["LN0438_MAAFHY1_R3"].pl.autocrop().pl.show(render_image=False, render_neighborhoods=True)

#plt.savefig("/g/saka/Tatjana/analysis/codex/scripts/testfig_nbh_5.png", dpi = 300, bbox_inches = "tight")
#plt.close()





