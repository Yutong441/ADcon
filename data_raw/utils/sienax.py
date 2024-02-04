# extract volume info from sienax report
import os
import re
import shutil
import subprocess
import nibabel as nib
from fsl.data.image import Image
import pandas as pd


def find_vol_per_line(string, feature):
    if re.search("^"+feature, string):
        new_str = re.sub(feature, "", string)
        new_str = re.sub("\\(.*\\)", "", new_str)
        new_str = re.sub(" ", ",", new_str)
        new_str = re.sub("\n", "", new_str)
        new_str = new_str.split(",")
        new_str = [i for i in new_str if i != ""]
        return new_str
    else:
        return None


def cal_csf_vol(sienax_dir):
    csf_seg = Image(sienax_dir+"/I_stdmaskbrain_pve_0.nii.gz")
    voxel_dim = csf_seg.header["pixdim"]
    voxel_vol = voxel_dim[1]*voxel_dim[2]*voxel_dim[3]
    csf_vol = voxel_vol*csf_seg.data.sum()
    return csf_vol


def extract_vol(filepath):
    sienax_dir = os.path.dirname(filepath)
    with open(filepath, "r") as f:
        text = f.readlines()

    vol = {"GREY": [], "WHITE": [], "BRAIN": [], "pgrey": [], "vcsf": []}
    for i in text:
        if re.match("^VSCALING.*", i):
            scaling = re.sub("^VSCALING", "", i)
            scaling = float(re.sub(" ", "", scaling))

        for j in vol:
            out = find_vol_per_line(i, j)
            if out is not None:
                vol[j] = [float(k) for k in out]

    csf_vol = cal_csf_vol(sienax_dir)
    vol["CSF"] = [csf_vol*scaling, csf_vol]
    vol["TOTAL"] = [vol["GREY"][0] + vol["WHITE"][0] + vol["CSF"][0],
                    vol["GREY"][1] + vol["WHITE"][1] + vol["CSF"][1]]

    vol = pd.DataFrame(vol)
    if vol.shape[0] == 2:
        vol.index = ["normalized", "unnormalized"]

    vol.to_csv(sienax_dir+"/volumes.csv")


def img2nii(filepath, save_path, header_path=None):
    '''Convert .img.gz file to nifti'''
    save_dir = os.path.dirname(save_path)
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    new_path = save_dir+"/"+os.path.basename(filepath)

    if re.match(".*.img.gz$", filepath):
        if header_path is None:
            header_path = re.sub(".img.gz$", ".hdr", filepath)
        new_header = re.sub(".img.gz$", ".hdr", new_path)
        decomp = True
    else:
        if header_path is None:
            header_path = re.sub(".img$", ".hdr", filepath)
        new_header = re.sub(".img$", ".hdr", new_path)
        decomp = False

    shutil.copyfile(filepath, new_path)
    shutil.copyfile(header_path, new_header)

    if decomp:
        cmd = ["gunzip", new_path]
        subprocess.run(cmd, stderr=subprocess.PIPE)

    img_nib = nib.load(re.sub(".gz$", "", new_path))
    new_img = nib.Nifti1Image(img_nib.get_fdata(), affine=img_nib.affine,
                              header=img_nib.header)
    new_img.to_filename(save_path)

    os.remove(re.sub(".gz$", "", new_path))
    os.remove(new_header)


def sienax_all(infile):
    sienax_dir = re.sub(".nii.gz$", "_sienax", infile)
    if not os.path.exists(sienax_dir+"/report.html"):
        cmd = ["sienax", infile, "-r", "-d"]
        subprocess.run(cmd, stderr=subprocess.PIPE)
        extract_vol(sienax_dir+"/report.sienax")
