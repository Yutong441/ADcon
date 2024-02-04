# plot p value along each tract in a template
import re
import os
import numpy as np
import nibabel as nib
import pandas as pd


def plot_tract(tract_val, tract_mask_dir, tract_conv_path,
               save_img_path, display_col="padj"):
    tract_conv = pd.read_csv(tract_conv_path)

    one_nib = nib.load(tract_mask_dir+"/"+os.listdir(tract_mask_dir)[0])
    affine = one_nib.affine
    all_tracts = np.zeros(one_nib.shape)
    for i in os.listdir(tract_mask_dir):
        ID = re.sub(".nii.gz", "", i)
        one_tract = nib.load(tract_mask_dir+"/"+i).get_fdata()
        one_label = tract_conv[tract_conv["abbrev"] == ID]
        one_label = one_label["full_name"].values[0]
        one_label = re.sub("_", " ", one_label)
        val_col = tract_val[tract_val["Module"] == one_label]
        if len(val_col) > 0:
            val_col = val_col[display_col].values[0]
            if val_col < 0.05:
                if display_col == "p" or display_col == "padj":
                    val_col = -np.log(val_col + 0.001)/np.log(10)
                if np.isnan(val_col):
                    val_col = 0
                all_tracts += one_tract*val_col

    affine = np.eye(4)
    affine[0, 0] = -1
    affine[1, 1] = -1
    all_tracts = nib.Nifti1Image(all_tracts, affine=affine)
    all_tracts.to_filename(save_img_path)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--df_path', type=str,
                        default="/home/yutong/Downloads/tracts_AD.csv")
    parser.add_argument('--save_dir', type=str,
                        default="/home/yutong/Data/ADcon/VLSM/tracts/")
    tract_mask_dir = "/home/yutong/Data/tract_data/tract_mask"
    tract_conv_path = "/home/yutong/Data/tract_data/fiber_abbrev.csv"
    args = parser.parse_args()

    tract_val = pd.read_csv(args.df_path)
    all_groups = pd.unique(tract_val["Group"])
    for i in all_groups:
        one_df = tract_val[tract_val["Group"] == i]
        plot_tract(one_df, tract_mask_dir, tract_conv_path,
                   args.save_dir+"/tract_"+i+".nii.gz",
                   display_col="padj")
