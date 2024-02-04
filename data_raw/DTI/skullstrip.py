import re
import os
import glob
import shutil
import subprocess
import pkg_resources

import numpy as np
import ants
import nibabel as nib
from fsl.wrappers import bet


def bash_in_python(bashcommand):
    ID = ''.join([str(i) for i in np.random.choice(9, 10)])
    sh_script = 'tmp_test'+ID+'.sh'
    with open(sh_script, 'w') as f:
        f.write('#!/bin/bash \n')
        f.write(bashcommand)
    subprocess.run(['chmod', '+x', sh_script], stderr=subprocess.PIPE)
    subprocess.call('./'+sh_script)
    os.remove(sh_script)


def degibbs(img, tmp_root="./tmp", nthreads=4):
    if not os.path.exists(tmp_root):
        os.mkdir(tmp_root)
    ID = ''.join([str(i) for i in np.random.choice(9, 20)])
    tmp_dir = tmp_root+"/"+ID
    os.mkdir(tmp_dir)
    ants.image_write(img, tmp_dir+"/original.nii.gz")
    cmd = ["mrdegibbs", "-nthreads", str(nthreads), "-force",
           tmp_dir+"/original.nii.gz", tmp_dir+"/degibbs.nii.gz"]
    subprocess.run(cmd, stderr=subprocess.PIPE)
    out = ants.image_read(tmp_dir+"/degibbs.nii.gz")
    shutil.rmtree(tmp_dir)
    return out


def ants_mask(img, mask):
    image_out = img.clone() * 0
    image_out[mask == 1] = img[mask == 1]
    return image_out


def skullstrip_reg(input_path, output_path, MNI=None):
    '''UKB technique: warp T1 to MNI non-linearly'''
    if MNI is None:
        MNI = pkg_resources.resource_filename(
                __name__, "../data/template/MNI.nii.gz")
    MNI_head = re.sub(".nii.gz", "_head.nii.gz", MNI)
    MNI_mask = re.sub(".nii.gz", "_mask.nii.gz", MNI)
    save_dir = os.path.dirname(output_path)
    template = ants.image_read(MNI_head)
    mask = ants.image_read(MNI_mask)
    img = ants.image_read(input_path)
    tx = ants.registration(template, img,
                           type_of_transform="SyN",
                           outprefix=save_dir+"/syn")
    mask_reg = ants.apply_transforms(img, mask,
                                     tx["invtransforms"])
    img_ss = ants_mask(img, mask_reg)
    for i in glob.glob(save_dir+"/syn*"):
        os.remove(i)
    ants.image_write(img_ss, output_path)
    ants.image_write(mask_reg, re.sub(".nii.gz", "_mask.nii.gz",
                                      output_path))


def skullstrip(input_path, output_path, method="hd-bet", MNI=None):
    '''
    Choose one of the following methods:
        `hd-bet`: deep learning tool for skullstripping of a variety of MRI
        modalities, but takes a long time on cpu
        `bet`: FSL BET
        `robex`
        `optibet`: use FNIRT to improve BET accuracy
        `UKB`: in UKB, T1 head is warped non-linearly to MNI head, then brain
        mask in MNI is inverse warped to T1 native space.
        `none`: no skullstripping is performed
    '''
    if not os.path.exists(output_path):
        if method == "hd-bet":
            cmd = ["hd-bet", "--input", input_path,
                   "-mode", "fast", "-device", "cpu", "-tta", "0",
                   "-o", output_path]
            bash_in_python(" ".join(cmd))
        elif method == "bet":
            bet(input_path, output_path, mask=True)
        elif method == "bet_dwi":
            bet(input_path, output_path, mask=True)
            cmd = ["bet", input_path, output_path, "-f", "0.3", "-R", "-m"]
            bash_in_python(" ".join(cmd))
        elif method == "robex":
            from pyrobex.robex import robex
            img = nib.load(input_path)
            stripped, mask = robex(img)
            stripped.to_filename(output_path)
            mask.to_filename(re.sub(".nii.gz", "_mask.nii.gz", output_path))
        elif method == "optibet":
            script_dir = os.path.dirname(__name__)
            cmd = ["bash", script_dir+"/optiBET.sh", "-i", input_path,
                   "-o", output_path]
            bash_in_python(" ".join(cmd))
        elif method == "UKB":
            skullstrip_reg(input_path, output_path, MNI)
        elif method == "none":
            shutil.copyfile(input_path, output_path)
            img_nib = nib.load(input_path)
            mask = np.ones(img_nib.shape)
            mask_nib = nib.Nifti1Image(mask, affine=img_nib.affine,
                                       header=img_nib.header)
            mask_nib.to_filename(re.sub(".nii.gz", "_mask.nii.gz",
                                        output_path))


def skullstrip_mdwi(input_path, output_path, b0_path, method="bet"):
    save_dir = os.path.dirname(output_path)
    if save_dir == "":
        save_dir = "."
    corr = ants.image_read(input_path)
    one_corr = ants.slice_image(corr, axis=3, idx=0)
    mdwi_ants = ants.from_numpy(corr.mean(3), origin=one_corr.origin,
                                spacing=one_corr.spacing,
                                direction=one_corr.direction)
    del corr, one_corr
    ants.image_write(mdwi_ants, save_dir+"/mdwi.nii.gz")
    del mdwi_ants
    skullstrip(save_dir+"/mdwi.nii.gz", output_path, method=method)
    b0_img = ants.image_read(b0_path)
    mask = ants.image_read(re.sub(".nii.gz", "_mask.nii.gz", output_path))
    b0_ss = ants.mask_image(b0_img, mask)
    ants.image_write(b0_ss, output_path)
    os.remove(save_dir+"/mdwi.nii.gz")
