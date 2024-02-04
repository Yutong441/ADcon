import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
from dipy.io.image import load_nifti_data
from dipy.io.streamline import load_tractogram
from dipy.tracking.streamline import select_random_set_of_streamlines
from dipy.tracking.utils import density_map
from dipy.tracking.streamline import length


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


def get_density(tract_path):
    stream, dimensions, _ = get_streamline(tract_path)
    stream_vox = select_random_set_of_streamlines(stream, 1000)
    density = density_map(stream_vox, np.eye(4), dimensions)
    return density


def to_cubic(img):
    ''' Pad a 3D image with zeros so that all 3 dimensions are equal '''
    img_shape = img.shape
    max_dim = np.max(np.array(img_shape))
    out = img.copy()
    for index, i in enumerate(img_shape):
        padding_left = list(out.shape)
        padding_left[index] = (max_dim - i)//2
        padding_right = list(out.shape)
        padding_right[index] = (max_dim - i)//2 + (max_dim - i) % 2
        out = np.concatenate([np.zeros(padding_left), out,
                              np.zeros(padding_right)], axis=index)
    return out


def plot3view(img, template, save_path=None):
    '''
    Plot axial, sagittal, and coronal MIP of `img` on the mid-slice of
    `template`
    '''
    assert len(img.shape) == 3
    assert len(template.shape) == 3
    H, W, D = template.shape

    plt.subplot(1, 3, 1)
    plt.imshow(template[..., D//2], cmap='gray')
    plt.imshow(img.max(2), alpha=0.5)
    plt.gca().set_title('axial')
    plt.axis('off')

    plt.subplot(1, 3, 2)
    plt.imshow(template[:, W//2, :], cmap='gray')
    plt.imshow(img.max(1), alpha=0.5)
    plt.gca().set_title('sagittal')
    plt.axis('off')

    plt.subplot(1, 3, 3)
    plt.imshow(template[H//2], cmap='gray')
    plt.imshow(img.max(0), alpha=0.5)
    plt.gca().set_title('coronal')
    plt.axis('off')
    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        plt.close()


def show_tracts(tract_path, save_path=None,
                MNI_path='tract_data/atlas/MNI.nii.gz'):
    density = get_density(tract_path)
    t1_data = load_nifti_data(MNI_path)
    plot3view(density, t1_data, save_path)
