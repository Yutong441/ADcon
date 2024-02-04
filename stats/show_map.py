import re
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import nibabel as nib
import pandas as pd


def disp_overlay(img, overlay, mask=None, save_path=None, width=6,
                 height=2, magni=6, exclude_slices=40, reverse_color=False,
                 minval=None, maxval=None):
    N = width*height
    D = img.shape[2]
    seq = np.linspace(exclude_slices, D - exclude_slices, num=N)
    seq = [int(i) for i in seq]

    if not reverse_color:
        cmap = mpl.cm.get_cmap("rainbow").copy()
    else:
        cmap = mpl.cm.get_cmap("rainbow_r").copy()
    cmap.set_under(color='black')
    mpl.rcParams.update({'font.size': 22})

    if mask is None:
        # mask = np.ones(img.shape)
        mask = (overlay != 0).astype(float)

    fig, axes = plt.subplots(height, width,
                             figsize=(magni*width, magni*height))
    if maxval is None:
        maxval = (mask*overlay).max()
    minval = (mask*overlay).min() + 0.001

    overlay *= mask
    overlay += (1 - mask)*(minval - 0.001)
    for index, i in enumerate(seq):
        axes[index // width, index % width].imshow(
                np.rot90(img[::-1, :, i], axes=(1, 0)), cmap="gray")
        im = axes[index // width, index % width].imshow(
            np.rot90(overlay[::-1, :, i], axes=(1, 0)), cmap=cmap, alpha=0.5,
            vmin=minval, vmax=maxval)
        axes[index // width, index % width].axis("off")
        axes[index // width, index % width].set_title("z={}".format(int(i)))

    plt.colorbar(im, ax=axes.ravel().tolist())
    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


# --------------------voxel wise Cox regression--------------------
template_path = "/home/yutong/Data/tract_data/template/MNI_header.nii.gz"
for i in ["EOAD", "LOAD"]:
    # for j in ["AD", "EF", "LAN", "MEM", "VS"]:
    for j in ["A5"]:
        var = i+"_"+j
        pval_path = "/home/yutong/Data/ADcon/VLSM/{}/pval_img.nii.gz".format(var)
        save_path_pval = "/home/yutong/Data/ADcon/figures/{}_pval.jpg".format(var)

        MNI = nib.load(template_path).get_fdata()
        pval_img = nib.load(pval_path).get_fdata()
        disp_overlay(MNI, pval_img, mask=pval_img < 0.05,
                     save_path=save_path_pval, reverse_color=True, minval=0.001)

var = "EOLO"
pval_path = "/home/yutong/Data/ADcon/VLSM/{}/pval_img.nii.gz".format(var)
save_path_pval = "/home/yutong/Data/ADcon/figures/{}_pval.jpg".format(var)

MNI = nib.load(template_path).get_fdata()
pval_img = nib.load(pval_path).get_fdata()
disp_overlay(MNI, pval_img, mask=pval_img < 0.05,
             save_path=save_path_pval, reverse_color=True, minval=0.001)

# --------------------frequency map--------------------
img_path = "/home/yutong/Data/ADcon/VLSM/EOAD_AD/freq_map.nii.gz"
save_path = "/home/yutong/Data/ADcon/figures/frequency_EOAD.jpg"
MNI = nib.load(template_path).get_fdata()
img = nib.load(img_path).get_fdata()
disp_overlay(MNI, img/120, save_path=save_path)

types = ["EOAD_2D", "EOAD_3D", "LOAD_2D", "LOAD_3D"]
root = "/home/yutong/Data/ADcon/"
for i in types:
    img_path = root+"/VLSM/{}/freq_map.nii.gz".format(i)
    save_path = root+"/figures/frequency/{}.jpg".format(i)
    df_path = root+"/VLSM/labels/{}.csv".format(i)
    N = pd.read_csv(df_path).shape[0]
    MNI = nib.load(template_path).get_fdata()
    img = nib.load(img_path).get_fdata()
    disp_overlay(MNI, img/N, save_path=save_path, maxval=0.9)

# --------------------tract map--------------------
img_dir = "/home/yutong/Data/ADcon/VLSM/tracts/"
save_dir = "/home/yutong/Data/ADcon/figures/tracts"
MNI = nib.load(template_path).get_fdata()

for i in os.listdir(img_dir):
    img = nib.load(img_dir+"/"+i).get_fdata()
    if img.sum() > 0:
        print(i)
        disp_overlay(MNI, img,
                     save_path=save_dir+"/"+re.sub(".nii.gz", ".jpg", i))
