# obtain connectivity matrix in normal and abnormal brain scans
import re
import os
import numpy as np
import pandas as pd
import nibabel as nib
from dipy.io.image import load_nifti_data
from dipy.tracking import utils
from dipy.io.streamline import load_tractogram
from dipy.tracking.streamline import length
from cdipy import conn_mat as CCM


def remove_bound(x, dimen, axis=2):
    keep = (x[:, axis] < (dimen[axis] - 0.5)) * ((x[:, axis] > 0))
    return x[keep]


def get_streamline(tract_path, reference=None, check=False,
                   upper_thres=300, lower_thres=30):
    '''
    Obtain streamlins from trk files
    Reference:
    https://dipy.org/documentation/1.5.0/examples_built/streamline_formats/#example-streamline-formats
    '''
    if reference is None:
        reference = 'same'
    else:
        reference = nib.load(reference)

    tract_ob = load_tractogram(tract_path, reference, bbox_valid_check=check)
    tract_ob.remove_invalid_streamlines()
    # Convert to world coordinates to obtain streamline length
    tract_ob.to_voxmm()
    tract_ob.to_corner()
    stream = tract_ob.streamlines
    affine, dimensions, voxel_sizes, voxel_order = tract_ob.space_attributes
    len_stream = np.array(length(stream))
    keep_lines = (len_stream > lower_thres) * (len_stream < upper_thres)

    # Switch to voxel coordinates to match to atlas
    tract_ob.to_vox()
    stream = tract_ob.streamlines[keep_lines]
    # remove the points that nearly touch the border
    stream = [np.array(remove_bound(i, dimensions),
                       dtype=np.float64, order="C") for i in stream]
    dimensions = np.array(dimensions, dtype=np.uint16)
    voxel_sizes = np.array(voxel_sizes, dtype=np.float64)
    return stream, dimensions, voxel_sizes


def get_conn_mat(stream, vox, mask_path, lesion_path=None):
    '''
    Obtain connectivity matrix given a tractogram, parcellation atlas, and
    lesion mask
    Args:
        `stream`: list of streamline coordinates
        `vox`: voxel dimension
        `mask_path`: path to parcellation atlas
        `lesion_path`: path to lesion mask
    Return:
        M: N x N numpy ndarray
        nM: N x N numpy ndarray
        number of streamlines
    '''
    affine = np.eye(4, dtype=np.float64, order="C")

    if lesion_path is not None:
        lesion = load_nifti_data(lesion_path)
        stream = list(utils.target(stream, affine, lesion))
        del lesion

    labels = load_nifti_data(mask_path)
    labels = np.array(labels, dtype=np.uint16, order="C")
    M, nM = CCM.connectivity_matrix(stream, affine, labels,
                                    vox, labels.max(), 0)
    return np.array(M)[1:, 1:], np.array(nM)[1:, 1:], len(stream)


def conn_mat_all(tract_path, mask_path, lesion_path=None):
    all_M, all_nM = [], []
    all_tracks = sorted(os.listdir(tract_path))
    all_T = pd.DataFrame(np.zeros([len(all_tracks), 1]), index=all_tracks,
                         columns=["num"])
    for i in all_tracks:
        try:
            stream, dimensions, vox = get_streamline(tract_path+'/'+i)
            M, nM, T = get_conn_mat(stream, vox, mask_path, lesion_path)
            all_T.loc[i] = T
            all_M.append(M)
            all_nM.append(nM)
        except Exception:
            print("cannot compute connectivity matrix for tract " + i)

    all_M = np.stack(all_M, axis=-1)
    all_nM = np.stack(all_nM, axis=-1)
    # ignore the first row and column, which considers background pixels
    final_M = all_M.sum(axis=2)
    final_nM = all_nM.sum(axis=2)

    # convert to dataframe
    parcel = pd.read_csv(re.sub(".nii.gz$", ".csv", mask_path),
                         index_col=[0]).values[:, 0]
    final_M = pd.DataFrame(final_M, index=parcel, columns=parcel)
    final_nM = pd.DataFrame(final_nM, index=parcel, columns=parcel)
    return final_M, final_nM, all_T


def get_conn_mat_all(tract_path, mask_path, lesion_path=None, ref_path=None):
    """
    Obtain connectivity matrix given a lesion mask.
    If lesion mask is not supplied, the matrix of normal brain will be
    returned.
    Args:
        `tract_path`: directory to list of white matter tracts in .trk format
        `mask_path`: parcellation atlas
        `lesion_path`: lesion mask
        `ref_path`: reference directory
    Returns:
        `final_M`: connectivity matrix
        `final_nM`: connectivity matrix normalized by streamline length
        `tract_num`: number of streamlines in each white matter tract
    """
    mask_ID = re.sub(".nii.gz$", "", os.path.basename(mask_path))
    if ref_path is not None:
        final_M_path = ref_path+"/ref_"+mask_ID+"_count.csv"
        final_nM_path = ref_path+"/ref_"+mask_ID+"_ncount.csv"
        tract_num_path = ref_path+"/ref_num.csv"
        ref_file = os.path.exists(final_M_path)
        if not os.path.exists(ref_path):
            os.mkdir(ref_path)
    else:
        ref_file = False

    if not ref_file or lesion_path is not None:
        final_M, final_nM, tract_num = conn_mat_all(tract_path, mask_path,
                                                    lesion_path)
        # if not ref_file:
        #     final_M.to_csv(final_M_path)
        #     final_nM.to_csv(final_nM_path)
        #     tract_num.to_csv(tract_num_path)

    elif ref_file and lesion_path is None:
        final_M = pd.read_csv(final_M_path, index_col=[0])
        final_nM = pd.read_csv(final_nM_path, index_col=[0])
        tract_num = pd.read_csv(tract_num_path, index_col=[0])
    return final_M, final_nM, tract_num


def get_tract_num(stream, vox, dimen, tract_dir, threshold=0.8):
    '''
    Obtain the number of streamlines in each white matter tract
    Args:
        `stream`: list of streamline coordinates
        `vox`: voxel dimension
        `tract_dir`: directory of storing the mask of every white matter tract
        registered to the native space
        `threshold`: a streamline is inside a white matter tract if how much
        percent of it is covered.
    '''
    affine = np.eye(4, dtype=np.float64, order="C")
    all_tracts = sorted(os.listdir(tract_dir))
    tract_stats = np.zeros([len(all_tracts), 4])

    for index, i in enumerate(all_tracts):
        tract_mask = load_nifti_data(tract_dir+"/"+i)
        tract_mask = np.array(tract_mask, dtype=np.uint16, order="C")
        num, fiber_len = CCM.target(stream, affine, tract_mask, vox,
                                    threshold)
        ROI = tract_mask.sum()
        dense = num/ROI
        dense_len = fiber_len/ROI

        tract_stats[index, 0] = num
        tract_stats[index, 1] = fiber_len
        tract_stats[index, 2] = dense
        tract_stats[index, 3] = dense_len

    features = ["num", "total length", "density", "normalized length"]
    tract_stats = pd.DataFrame(tract_stats, index=all_tracts,
                               columns=features)
    img = CCM.fiber_density_map(stream, affine, dimen)
    return tract_stats, img
