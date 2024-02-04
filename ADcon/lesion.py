import os
import shutil
import numpy as np
import ants
import SimpleITK as sitk
import pandas as pd
import conn_metric as ME


def read_ants(img_path):
    img = sitk.ReadImage(img_path)
    direction = np.array(img.GetDirection()).reshape(3, 3)
    spacing = img.GetSpacing()
    origin = img.GetOrigin()
    np_img = sitk.GetArrayFromImage(img).astype(float)
    img_ant = ants.from_numpy(np.transpose(np_img, (2, 1, 0)),
                              origin=origin,
                              spacing=spacing, direction=direction)
    return img_ant


def syn_trans(wdir, MNI, tmp_dir):
    T1 = ants.image_read(wdir+"/DL_yc/T1_MNI.nii.gz")
    tx = ants.registration(MNI, T1, type_of_transform="SyN",
                           outprefix=tmp_dir+"/ants", random_seed=100)
    return tx["fwdtransforms"]


def reg_WMH(wdir, save_path, wm_template, mni_template):
    if os.path.exists(wdir+"/DL_yc/tissue_MNI.nii.gz") and os.path.exists(
            wdir+"/DL_yc/WMH_MNI.nii.gz"):
        MNI = read_ants(mni_template)
        ID = "".join([str(i) for i in np.random.choice(10, 10)])
        tmp_dir = os.path.dirname(save_path)+"/tmp"+ID+"/"

        # compute SyN transform
        os.mkdir(tmp_dir)
        transform_files = syn_trans(wdir, MNI, tmp_dir)

        # warp lesion mask to MNI space
        # select correct mask
        if os.path.exists(wdir+"/DL_yc/WMH_MNIH.nii.gz"):
            img_path = wdir+"/DL_yc/WMH_MNIH.nii.gz"
        else:
            img_path = wdir+"/DL_yc/WMH_MNI.nii.gz"
        lesion_mask = read_ants(img_path)
        lesion_mask = ants.apply_transforms(
            MNI, lesion_mask, transform_files, interpolator="genericLabel")

        # warp WM mask to MNI space
        tissue_mask = ants.image_read(wdir+"/DL_yc/tissue_MNI.nii.gz")
        tissue_mask = ants.threshold_image(tissue_mask, low_thresh=3,
                                           binary=True)
        tissue_mask = ants.apply_transforms(
            MNI, tissue_mask, transform_files, interpolator="genericLabel")

        tissue_MNI = ants.image_read(wm_template)

        prefix = tmp_dir+"/mask"
        tx = ants.registration(tissue_MNI, tissue_mask,
                               type_of_transform="SyN",
                               outprefix=prefix, random_seed=100)
        lesion_MNI = ants.apply_transforms(
            tissue_MNI, lesion_mask, tx["fwdtransforms"],
            interpolator="genericLabel")

        shutil.rmtree(tmp_dir)
        ants.image_write(lesion_MNI, save_path)


def lesion_metric(final_dir, atlas_dir):
    all_me, all_nodes = [], []
    if not os.path.exists(final_dir+"/metric"):
        os.mkdir(final_dir+"/metric")
    if os.path.exists(final_dir+"/WMH_reg.nii.gz"):
        metrics, nodes = ME.get_all_metrics(
            final_dir+"/WMH_reg.nii.gz", final_dir+"/metric/", atlas_dir)
        all_me.append(metrics)
        all_nodes.append(nodes)

    all_me = pd.concat(all_me, axis=0)
    all_nodes = pd.concat(all_nodes, axis=0)
    all_me.to_csv(final_dir+"/graph_metrics.csv")
    all_nodes.to_csv(final_dir+"/node_metrics.csv")


def lesion_disconn(root_dir, save_dir, atlas_dir, num, arrayID):
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    for index, i in enumerate(os.listdir(root_dir)[0:1]):
        if index % num == arrayID:
            if not os.path.exists(save_dir+"/"+i):
                os.mkdir(save_dir+"/"+i)

            if not os.path.exists(save_dir+"/"+i+"/WMH_reg.nii.gz"):
                reg_WMH(root_dir+"/"+i, save_dir+"/"+i+"/WMH_reg.nii.gz",
                        atlas_dir+"/template/MNI_WM.nii.gz",
                        atlas_dir+"/template/MNI.nii.gz")
            if not os.path.exists(save_dir+"/"+i+"/metric/tract_num.csv"):
                if os.path.exists(save_dir+"/"+i+"/WMH_reg.nii.gz"):
                    lesion_metric(save_dir+"/"+i, atlas_dir)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--root_dir', type=str)
    parser.add_argument('--save_dir', type=str)
    parser.add_argument('--atlas_dir', type=str)
    parser.add_argument("--num", type=int, default=1)
    parser.add_argument("--arrayID", type=int, default=0)
    args = parser.parse_args()
    lesion_disconn(args.root_dir, args.save_dir, args.atlas_dir, args.num,
                   args.arrayID)
