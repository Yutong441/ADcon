# white matter hyperintensity segmentation using Hypermapp3r
import os
import shutil
import subprocess
import numpy as np


def bash_in_python(bashcommand):
    ID = ''.join([str(i) for i in np.random.choice(9, 10)])
    sh_script = 'tmp_test'+ID+'.sh'
    with open(sh_script, 'w') as f:
        f.write('#!/bin/bash \n')
        f.write(bashcommand)
    subprocess.run(['chmod', '+x', sh_script], stderr=subprocess.PIPE)
    subprocess.call('./'+sh_script)
    os.remove(sh_script)


def hypermap(img_dir, tmp_root, env_path):
    if not os.path.exists(img_dir+"/WMH_MNI.nii.gz"):
        if not os.path.exists(tmp_root):
            os.mkdir(tmp_root)
        ID = ''.join([str(i) for i in np.random.choice(9, 10)])
        tmp_dir = tmp_root+"/"+ID
        os.mkdir(tmp_dir)
        for i in ["T1", "FLAIR", "mask"]:
            shutil.copyfile(img_dir+"/{}_MNI.nii.gz".format(i),
                            tmp_dir+"/{}_MNI.nii.gz".format(i))

        ori_dir = os.getcwd()
        os.chdir(tmp_dir)
        bash_in_python(
            "source {}/venv_dwi/bin/activate;".format(env_path) +
            "hypermapper seg_wmh -fl FLAIR_MNI.nii.gz " +
            "-t1 T1_MNI.nii.gz -m mask_MNI.nii.gz -o WMH_MNI.nii.gz;" +
            "deactivate"
        )
        os.chdir(ori_dir)
        shutil.copyfile(tmp_dir+"/WMH_MNI.nii.gz", img_dir+"/WMH_MNI.nii.gz")
        shutil.rmtree(tmp_dir)

# def hypermap(img_dir, tmp_root, env_path, to_mni, mni, save_dir=None):
#     if save_dir is None:
#         save_dir = img_dir+"/MNI"
# 
#     if not os.path.exists(save_dir+"/WMH_MNI.nii.gz"):
#         if not os.path.exists(tmp_root):
#             os.mkdir(tmp_root)
#         ID = ''.join([str(i) for i in np.random.choice(9, 10)])
#         tmp_dir = tmp_root+"/"+ID
#         os.mkdir(tmp_dir)
# 
#         # register flair to t1
#         shutil.copyfile(img_dir+"/T1_reg.nii.gz", tmp_dir+"/T1.nii.gz")
#         shutil.copyfile(img_dir+"/FLAIR_reg.nii.gz", tmp_dir+"/FLAIR.nii.gz")
#         shutil.copyfile(img_dir+"/T1_ss_mask.nii.gz", tmp_dir+"/mask.nii.gz")
# 
#         ori_dir = os.getcwd()
#         os.chdir(tmp_dir)
#         bash_in_python(
#             "source {}/venv_dwi/bin/activate;".format(env_path) +
#             "hypermapper seg_wmh -fl FLAIR.nii.gz " +
#             "-t1 T1.nii.gz -m mask.nii.gz -o WMH.nii.gz;" +
#             "deactivate"
#         )
#         os.chdir(ori_dir)
# 
#         wmh = ants.image_read(tmp_dir+"/WMH.nii.gz")
#         wmh = ants.apply_transforms(mni, wmh, [to_mni])
#         ants.image_write(wmh, tmp_dir+"/WMH_MNI.nii.gz")
#         shutil.copyfile(tmp_dir+"/WMH_MNI.nii.gz", save_dir+"/WMH_MNI.nii.gz")
#         shutil.copyfile(tmp_dir+"/WMH.nii.gz", img_dir+"/WMH_T1.nii.gz")
#         shutil.rmtree(tmp_dir)


def all_segs(img_root, tmp_dir, num=1, arrayID=0):
    for index, i in enumerate(sorted(os.listdir(img_root))):
        if index % num == arrayID:
            if os.path.exists(img_root+"/"+i+"/FLAIR_MNI.nii.gz"):
                if not os.path.exists(img_root+"/"+i+"/WMH_MNI.nii.gz"):
                    hypermap(img_root+"/"+i, tmp_dir)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--img_dir", type=str)
    parser.add_argument("--tmp_dir", type=str)
    parser.add_argument("--num", type=int)
    parser.add_argument("--arrayID", type=int)
    args = parser.parse_args()
    all_segs(args.img_dir, args.tmp_dir, args.num, args.arrayID)
