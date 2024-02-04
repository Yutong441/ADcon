import nibabel as nib
import pandas as pd


def pipeline(wdir, save_path, df_path=None):
    if df_path is None:
        df_path = "/home/yutong/Data/postICH/cleaned_labels/all_sel.csv"

    all_df = pd.read_csv(df_path)
    all_IDs = all_df["anoID"].values
    # all_files = os.listdir(wdir)
    template = nib.load(wdir+"/"+all_IDs[0]+".nii.gz")
    templ_img = template.get_fdata().astype(float)

    for i in all_IDs:
        new_img = nib.load(wdir+"/"+i+".nii.gz").get_fdata().astype(float)
        templ_img += new_img

    out = nib.Nifti1Image(templ_img, affine=template.affine,
                          header=template.header)
    out.to_filename(save_path)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--wdir', type=str,
                        default="/home/yutong/Data/ICHmap/reg_lesion/")
    parser.add_argument('--save_path', type=str,
                        default="/home/yutong/Data/ICHmap/results/freq_maps.nii.gz")
    args = parser.parse_args()
    pipeline(args.wdir, args.save_path)
