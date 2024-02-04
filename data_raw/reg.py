import re
import os
import json
import glob
import shutil
import time
import pkg_resources

import numpy as np
import ants
import nibabel as nib

from data_raw.DTI import skullstrip as SS
from data_raw.utils.WMH_seg import hypermap
from data_raw.utils.sienax import sienax_all


def register_ss(img_paths, template_path, save_dir, N4=True):
    ''' Register all modalities except DWI to T1, and skullstrip T1 '''
    template = ants.image_read(template_path)
    template = SS.degibbs(template, tmp_root=save_dir)
    template = ants.n4_bias_field_correction(template)
    template.to_filename(save_dir+"/T1_reg.nii.gz")

    # hd-bet produces better results than BET or RODEX, but takes longer
    SS.skullstrip(save_dir+"/T1_reg.nii.gz", save_dir+"/T1_ss.nii.gz",
                  method="hd-bet")
    sienax_all(save_dir+"/T1_reg.nii.gz")

    # register other modalities to T1
    for key, img_path in img_paths.items():
        if os.path.exists(img_path):
            new_path = img_path
            img = ants.image_read(new_path)
            img = SS.degibbs(img, tmp_root=save_dir)
            if len(img.shape) == 4:
                img = ants.slice_image(img, axis=3, idx=0)
            if N4:
                img = ants.n4_bias_field_correction(img)

            ants.image_write(img, save_dir+"/"+key+"_nat.nii.gz")
            tx = ants.registration(template, img,
                                   type_of_transform="Rigid",
                                   outprefix=save_dir+"/"+key)

            ants.image_write(tx['warpedmovout'],
                             save_dir+"/"+key+"_reg.nii.gz")
            # rename and save transform
            os.rename(tx["fwdtransforms"][0],
                      save_dir+"/"+key+"_to_T1_rigid.mat")


def load_template(img_path):
    '''Crop images in MNI space into 160 x 192 x 128 images'''
    img = nib.load(img_path)
    if img.shape[0] == 182 and img.shape[1] == 218 and img.shape[2] == 182:
        img_new = nib.Nifti1Image(img.get_fdata()[12:172, 16:208, 22:150],
                                  affine=img.affine, header=img.header)
        return ants.from_nibabel(img_new)
    else:
        return ants.image_read(img_path)


def standardize_img(img, mask):
    masked_img = np.ma.masked_array(img.get_fdata(),
                                    mask=mask.get_fdata() == 0)
    masked_img = (masked_img - masked_img.mean())/masked_img.std()
    masked_img = masked_img.clip(-2, 2)
    masked_img = (masked_img + 2)/4
    stand_img = np.ma.filled(masked_img, fill_value=0)
    stand_img = nib.Nifti1Image(stand_img, header=img.header,
                                affine=img.affine)
    return stand_img


def warp_MNI(save_dir, others_dict, MNI_path):
    if not os.path.exists(save_dir+"/MNI"):
        os.mkdir(save_dir+"/MNI")

    # rigid registration of skullstripped T1 to MNI
    template = load_template(MNI_path)
    t1 = ants.image_read(save_dir+"/T1_ss.nii.gz")
    tx = ants.registration(template, t1, type_of_transform="Rigid",
                           outprefix=save_dir+"/MNI/warp")

    # save transforms
    os.rename(tx["fwdtransforms"][0], save_dir+"/MNI/T1_to_MNI_rigid.mat")

    # warp mask
    mask = ants.image_read(save_dir+"/T1_ss_mask.nii.gz",
                           pixeltype="unsigned int")
    mask_reg = ants.apply_transforms(
        template, mask, [save_dir+"/MNI/T1_to_MNI_rigid.mat"],
        interpolator="genericLabel")
    ants.image_write(mask_reg, save_dir+"/MNI/mask_MNI.nii.gz")
    mask_reg = ants.to_nibabel(mask_reg)

    # save standardized T1
    img_reg = standardize_img(ants.to_nibabel(tx["warpedmovout"]),
                              mask_reg)
    img_reg.to_filename(save_dir+"/MNI/T1_MNI.nii.gz")

    # register DWI
    coreg_mod = list(others_dict.keys())
    for i in coreg_mod:
        if os.path.exists(save_dir+"/"+i+"_reg.nii.gz"):
            img = ants.image_read(save_dir+"/"+i+"_reg.nii.gz")
            img = SS.ants_mask(img, mask)
            img_reg = ants.apply_transforms(
                template, img, save_dir+"/MNI/T1_to_MNI_rigid.mat")
            if i in coreg_mod:
                img_reg = standardize_img(ants.to_nibabel(img_reg),
                                          mask_reg)
            img_reg.to_filename(save_dir+"/MNI/"+i+"_MNI.nii.gz")

        # in direct registration, MNI registration took place before mask
        # registration. Thus, apply the brain mask
        elif os.path.exists(save_dir+"/MNI/"+i+"_MNI.nii.gz"):
            img = ants.image_read(save_dir+"/"+i+"_reg.nii.gz")
            img = SS.ants_mask(img, mask)
            ants.image_write(img, save_dir+"/MNI/"+i+"_MNI.nii.gz")


def analyses(save_dir, MNI_path, env_path):
    template = load_template(MNI_path)

    # tissue segmentation
    SS.bash_in_python("fast "+save_dir+"/T1_ss.nii.gz")
    os.rename(save_dir+"/T1_ss_seg.nii.gz",
              save_dir+"/tissue_T1.nii.gz")
    tissue_ss = ants.image_read(save_dir+"/tissue_T1.nii.gz")
    tissue_reg = ants.apply_transforms(
        template, tissue_ss, [save_dir+"/MNI/T1_to_MNI_rigid.mat"],
        interpolator="genericLabel")
    ants.image_write(tissue_reg, save_dir+"/MNI/tissue_MNI.nii.gz")

    # WMH segmentation
    if os.path.exists(save_dir+"/MNI/FLAIR_MNI.nii.gz"):
        hypermap(save_dir+"/MNI", save_dir+"/WMH_seg", env_path)


def pipeline_native(wdir, mask_path, save_dir):
    '''Inverse warp the T1 mask to native space of individual sequences'''
    mask = ants.image_read(mask_path)
    for i in glob.glob(wdir+"/*_nat.nii.gz"):
        key = re.sub("_nat.nii.gz$", "", os.path.basename(i))
        if key != "T1":
            img = ants.image_read(i)
            inv_mask = ants.apply_transforms(
                img, mask, [wdir+"/transform/"+key+"_to_T1_rigid.mat"],
                whichtoinvert=[True])
            masked = ants.mask_image(img, inv_mask)
            ants.image_write(masked, save_dir+"/"+key+"_native.nii.gz")


def mkdir(wdir):
    if not os.path.exists(wdir):
        os.mkdir(wdir)


def clean_dir(wdir):
    # store transform
    mkdir(wdir+"/transform")
    for i in glob.glob(wdir+"/MNI/T1_to_MNI*"):
        os.rename(i, wdir+"/transform/"+os.path.basename(i))
    for i in glob.glob(wdir+"/*_to_T1_rigid.mat"):
        shutil.copyfile(i, wdir+"/transform/"+os.path.basename(i))

    # store native space data
    mkdir(wdir+"/native")
    shutil.copyfile(wdir+"/tissue_T1.nii.gz",
                    wdir+"/native/tissue_native.nii.gz")
    shutil.copyfile(wdir+"/T1_ss_mask.nii.gz",
                    wdir+"/native/mask_native.nii.gz")
    shutil.copyfile(wdir+"/T1_ss.nii.gz", wdir+"/native/T1_native.nii.gz")
    pipeline_native(wdir, wdir+"/native/mask_native.nii.gz", wdir+"/native")

    # store ready-to-use image data
    if os.path.exists(wdir+"/DL_yc"):
        shutil.rmtree(wdir+"/DL_yc")
    os.rename(wdir+"/MNI", wdir+"/DL_yc")

    # stats
    mkdir(wdir+"/stats")
    shutil.copyfile(wdir+"/T1_reg_sienax/volumes.csv",
                    wdir+"/stats/brain_volumes.csv")

    # delete everything else
    for i in os.listdir(wdir):
        if not os.path.isdir(wdir+"/"+i):
            os.remove(wdir+"/"+i)

    shutil.rmtree(wdir+"/T1_reg_sienax")
    if os.path.exists(wdir+"/WMH_seg"):
        shutil.rmtree(wdir+"/WMH_seg")


def load_existing(save_dir):
    '''Load previous results from computationally heavy steps'''
    if not os.path.exists(save_dir+"/MNI"):
        os.mkdir(save_dir+"/MNI")
    if os.path.exists(save_dir+"/native/T1_native.nii.gz"):
        os.rename(save_dir+"/native/T1_native.nii.gz",
                  save_dir+"/T1_ss.nii.gz")
        os.rename(save_dir+"/native/mask_native.nii.gz",
                  save_dir+"/T1_ss_mask.nii.gz")
    if os.path.exists(save_dir+"/DL_yc/WMH_MNI.nii.gz"):
        os.rename(save_dir+"/DL_yc/WMH_MNI.nii.gz",
                  save_dir+"/MNI/WMH_MNI.nii.gz")
        shutil.rmtree(save_dir+"/DL_yc")


def pipeline(config, mni_path, others_dict,
             eddy="RUNDMC", dti_reg_method="indirect",
             dwi_ss_method="mdwi_bet", N4=True, env_path=None):
    start = time.time()
    proceed = [os.path.exists(config[i]) for i in ["T1"]]
    if sum(proceed) == 1:
        if not os.path.exists(config["save"]):
            os.mkdir(config["save"])
        load_existing(config["save"])
        register_ss(others_dict, config["T1"], config["save"], N4=N4)
        warp_MNI(config["save"], others_dict, mni_path)
        analyses(config["save"], mni_path, env_path)
        clean_dir(config["save"])
    else:
        print(config["T1"])
        print("T1, FLAIR or DWI do not exist, abort")
    print((time.time() - start)/61)


def modify_str(txt, ID):
    ''' translate regular expression in the data_regex.json file '''
    if "*8IDIDIDID" in txt:
        numer = re.sub("([a-z]|[A-Z])+", "", ID)
        ID_num = numer.zfill(8)
        return re.sub("8IDIDIDID", ID_num, txt)
    elif "IDIDIDID" in txt:
        return re.sub("IDIDIDID", ID, txt)


def _add_ID(img_list, root, ID):
    ilist = [root+"/"+modify_str(i, ID) for i in img_list]
    # account for regular expressions in json file
    for index, i in enumerate(ilist):
        if "*" in i:
            out = glob.glob(i)
            if len(out) > 0:
                ilist[index] = out[0]
    return ilist


def dwi_config(config: dict, root: str, ID: str, env_path: str):
    '''
    Convert the regex file pattern in config into specific files
    '''
    cf = {}
    cf["T1"] = _add_ID([config["T1"]], root, ID)[0]
    # cf["FLAIR"] = _add_ID([config["FLAIR"]], root, ID)[0]
    if "exception" in config.keys():
        cf["exception"] = config["exception"].split(",")
    else:
        cf["exception"] = []
    cf["save"] = root+"/"+config["save"]+"/"+ID
    return cf


class Reg_dataset:
    def __init__(self, **kwargs):
        self.opt = kwargs
        cf_name = pkg_resources.resource_filename(__name__,
                                                  "data_regex.json")
        with open(cf_name, "r") as f:
            cf = json.load(f)
        self.one_cf = cf[self.opt["data_name"]]
        if not os.path.exists(self.opt["root"]+"/"+self.one_cf["save"]):
            os.mkdir(self.opt["root"]+"/"+self.one_cf["save"])

    def run(self, num=1, arrayID=0):
        all_IDs = sorted(os.listdir(self.opt["root"]+"/"+self.one_cf["ID"]))
        for index, i in enumerate(all_IDs):
            # only run at a specific cpu (arrayID) out of a total of (num)
            # this is specified in a slurm script
            if index % num == arrayID:
                config = dwi_config(self.one_cf, self.opt["root"], i,
                                    self.opt["env_path"])
                if i not in config["exception"]:
                    if not os.path.exists(
                            config["save"]+"/DL_yc/T1_MNI.nii.gz"):
                        # for other non-DWI modalities co-registered with T1
                        others_dict = {}
                        for key, val in self.one_cf.items():
                            if key not in [*config.keys(), "ID"]:
                                others_dict[key] = _add_ID(
                                        [val], self.opt["root"], i)[0]

                        pipeline(config,
                                 mni_path=self.opt["mni_path"],
                                 others_dict=others_dict,
                                 eddy=self.opt["eddy"],
                                 dti_reg_method=self.opt["dti_reg_method"],
                                 dwi_ss_method=self.opt["dwi_ss_method"],
                                 env_path=self.opt["env_path"],
                                 N4=self.opt["N4"])


if __name__ == "__main__":
    MNI_path = os.environ["FSLDIR"] + \
        "/data/standard/MNI152_T1_1mm_brain.nii.gz"
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_name', type=str, default="SCANS_t0")
    parser.add_argument('--root', type=str)
    parser.add_argument("--mni", type=str, default=MNI_path)
    parser.add_argument("--eddy", type=str, default="RUNDMC")
    parser.add_argument("--dti_reg_method", type=str, default="indirect")
    parser.add_argument("--dwi_ss_method", type=str, default="bet")
    parser.add_argument("--no_N4", action="store_true")
    parser.add_argument("--env_path", type=str, default=".")
    parser.add_argument("--num", type=int, default=1)
    parser.add_argument("--arrayID", type=int, default=0)

    args = parser.parse_args()

    Reg = Reg_dataset(
            data_name=args.data_name,
            root=args.root,
            mni_path=args.mni,
            eddy=args.eddy,
            dti_reg_method=args.dti_reg_method,
            dwi_ss_method=args.dwi_ss_method,
            N4=not args.no_N4,
            env_path=args.env_path
    )
    Reg.run(args.num, args.arrayID)
