# quantify WMH and brain volumes
import os
import subprocess
import numpy as np
import pandas as pd
import nibabel as nib


def bash_in_python(cmd):
    out, err = unix_cmd(cmd)
    show_error(err)
    return out, err


def show_error(err):
    if len(err) > 0:
        print(err)


def unix_cmd(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
    out, err = p.communicate()
    return out, err


def custom_mask(arr, labels):
    out = np.zeros(arr.shape)
    for i in labels:
        out += arr == i
    return out > 0


def get_volumes(wdir, num_cpu=3):
    if not os.path.exists(wdir+"/segmentation/"):
        os.mkdir(wdir+"/segmentation/")

    if not os.path.exists(wdir+"/segmentation/seg.nii.gz"):
        bash_in_python(
            "mri_synthseg --i {}/native/T1_native.nii.gz".format(wdir) +
            " --o {}/segmentation/seg.nii.gz".format(wdir) +
            " --threads {} --fast".format(num_cpu)
        )

    out = {}

    seg_nib = nib.load(wdir+"/segmentation/seg.nii.gz")
    seg = seg_nib.get_fdata()
    pixdim = seg_nib.header["pixdim"]
    vox_dim = pixdim[1]*pixdim[2]*pixdim[3]
    # remove CSF
    seg[seg == 24] = 0
    WM_mask = custom_mask(seg, [2, 41, 7, 46, 16])
    out["TBV"] = (seg > 0).sum()*vox_dim
    out["WM_vol"] = WM_mask.sum()*vox_dim

    wmh_path = wdir+"/DL_yc/WMH_MNIH.nii.gz"
    if not os.path.exists(wmh_path):
        wmh_path = wdir+"/DL_yc/WMH_MNI.nii.gz"

    if os.path.exists(wmh_path):
        WMH_nib = nib.load(wmh_path)
        WMH = WMH_nib.get_fdata()
        res = WMH_nib.header["pixdim"][1:4]
        voxel_size = res[0]*res[1]*res[2]
        WMH_vol = WMH.sum()*voxel_size
        out["WMH_vol"] = WMH_vol

    out = pd.DataFrame.from_dict(out, columns=[os.path.basename(wdir)],
                                 orient="index")
    return out.T


def get_all_volumes(wdir, save_path, num=1, arrayID=0):
    all_out = []
    for index, i in enumerate(os.listdir(wdir)):
        if index % num == arrayID:
            if os.path.exists(wdir+"/"+i+"/stats/brain_volumes.csv"):
                all_out.append(get_volumes(wdir+"/"+i))

    all_out = pd.concat(all_out, axis=0)
    all_out.to_csv(save_path)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--wdir', type=str)
    parser.add_argument('--save_path', type=str)
    parser.add_argument("--num", type=int, default=1)
    parser.add_argument("--arrayID", type=int, default=0)
    args = parser.parse_args()
    get_all_volumes(args.wdir, args.save_path, args.num, args.arrayID)
