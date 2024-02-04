#cython: boundscheck=True
#cython: nonecheck=False
#cython: wraparound=False
#cython: infertypes=False
#cython: initializedcheck=False
import numpy as np
from libc.math cimport sqrt


cdef double c_length(double [:,::1] streamline, const double [:] weighting) nogil:
    '''
    This function is rewritten from `dipy.tracking.streamlinespeed.c_length` to
    calculate path length if a voxel is not isotropic.
    '''
    cdef:
        int i, j
        double out = 0.0
        double dn, sum_dn_sqr

    for i in range(1, streamline.shape[0]):
        sum_dn_sqr = 0.0
        for j in range(3):
            dn = streamline[i, j] - streamline[i-1, j]
            dn = dn * weighting[j]
            sum_dn_sqr += dn*dn

        out += sqrt(sum_dn_sqr)

    return out


def length(list streamlines, const double [:] weighting):
    cdef int L = len(streamlines)
    cdef int i
    cdef double one_l
    cdef double [:] streamline_length = np.zeros([L])
    cdef double [:,::1] one_stream

    for i in range(L):
        one_stream = streamlines[i]
        one_l = c_length(one_stream, weighting)
        streamline_length[i] = one_l
    return streamline_length


def get_ncount(dict group, int N, const double [:] weighting):
    '''
    Given a dictionary of streamlines, calculate the sum of inverse lengths of
    the streamlines connecting different regions.
    '''
    cdef int i, j
    cdef double sum_n
    cdef double [:,:] ncount = np.zeros([N, N])
    cdef double [:] streamline_length

    for i in range(N):
        for j in range(N):
            if (i, j) in group.keys():
                streamline_length = length(group[(i, j)], weighting)
                sum_n = 0
                for k in range(len(streamline_length)):
                    sum_n += 1/streamline_length[k]

                ncount[i, j] = sum_n 
                ncount[j, i] = sum_n
    return ncount


def _mapping_to_voxel(double [:,::1] affine):
    cdef double [:,::1] lin_T, inv_affine
    cdef double [:] offset
    cdef int i

    inv_affine = np.linalg.inv(affine)
    lin_T = inv_affine[:3, :3].T.copy()
    offset = inv_affine[:3, 3]
    for i in range(3):
        offset[i] += 0.5
    return lin_T, offset


cpdef _to_voxel_coordinates(double [:,::1] streamline, 
                          double [:,::1] lin_T, 
                          double [:] offset):
    ''' Dot product between streamline and lin_T plus offset '''
    cdef unsigned short [:,:] inds
    cdef int i, j, k
    cdef int npoint = len(streamline)
    cdef double sum_n

    inds = np.empty((npoint, 3), dtype=np.uint16)
    for i in range(npoint):
        for j in range(3):
            sum_n = 0
            for k in range(3):
                sum_n += streamline[i, k]*lin_T[k, j]

            inds[i, j] = <int>(sum_n + offset[j])
    return inds


def connectivity_matrix(list streamlines, double [:, ::1] affine, 
                        unsigned short [:, :, ::1] label_volume, 
                        double [:] weighting, int max_lab, int ends):
    '''
    This function is rewritten from `dipy.tracking.utils.connectivity_matrix`.
    It only works for symmetric matrix and considers the regions crossed by the
    entire streamline (`inclusive=True`).
    Furthermore, it does not return streamlines as output.

    Args:
        `streamlines`: list of streamline coordinates
        `affine`: voxel to physical space transform (typically identity matrix)
        `label_volume`: parcellation atlas
        `weighting`: voxel dimension
        `max_lab`: maximum label number

    Returns:
        `matrix`: count matrix
        `nmatrix`: sum of inverse streamline length
    '''
    cdef int sl, a, b, a2, b2, i, j, k, one_lab, label_len, line_len
    cdef double [:,::1] lin_T
    cdef double [:] offset
    cdef double [:,::1] one_line, new_one
    cdef unsigned short [:,::1] entire
    cdef unsigned short [::] labels, label_bool
    cdef list entireLabels

    cdef unsigned long [:,::1] matrix
    cdef double [:,::1] nmatrix
    cdef int mx = max_lab + 1

    matrix = np.zeros([mx, mx], dtype=np.uint64)
    nmatrix = np.zeros([mx, mx])
    labels = np.arange(mx, dtype=np.uint16)
    lin_T, offset = _mapping_to_voxel(affine)

    for sl in range(len(streamlines)):
        # Convert streamline to voxel coordinates
        one_line = streamlines[sl]
        line_len = len(one_line)
        if ends > 0:
            new_one = np.empty((ends*2, 3), dtype=np.float64)
            new_one[:ends, :] = one_line[:ends, :]
            new_one[ends:, :] = one_line[-ends:, :]
        else:
            new_one = one_line
        entire = _to_voxel_coordinates(new_one, lin_T, offset)

        # Create list of all labels streamline passes through
        label_bool = np.zeros([mx], dtype=np.uint16)
        for k in range(len(entire)):
            one_lab = label_volume[entire[k, 0], entire[k, 1], entire[k, 2]]
            label_bool[one_lab] = 1

        entireLabels = []
        for i in range(mx):
            if label_bool[i] == 1:
                entireLabels.append(labels[i])

        label_len = len(entireLabels)
        if label_len > 1:
            # Append all connection combinations with streamline number
            for a in range(label_len):
                for b in range(label_len):
                    if a > b:
                        a2 = entireLabels[a]
                        b2 = entireLabels[b]
                        matrix[a2, b2] += 1
                        nmatrix[a2, b2] += 1/c_length(one_line, weighting)

    # make a symmetric matrix
    with nogil:
        for i in range(mx):
            for j in range(mx):
                if i != j and matrix[i, j] == 0:
                    matrix[i, j] = matrix[j, i]
                    nmatrix[i, j] = nmatrix[j, i]

    return matrix, nmatrix



def target(list streamlines, double [:, ::1] affine, 
           unsigned short [:, :, ::1] target_mask, 
           double [:] weighting, double threshold):
    '''
    Calculate the number and combind length of the fibers passing through a
    lesion mask
    '''
    cdef int i, k, num_voxel, ROI, one_len
    cdef unsigned short one_label, state, thres
    cdef double [:,::1] one_str
    cdef unsigned short [:,::1] ind
    cdef double [:,::1] lin_T
    cdef double [:] offset
    # output
    cdef int num_fiber
    cdef double all_length

    lin_T, offset = _mapping_to_voxel(affine)
    num_fiber = 0
    all_length = 0.
    thres = 0
    for i in range(len(streamlines)):
        one_str = streamlines[i]
        ind = _to_voxel_coordinates(one_str, lin_T, offset)
        one_len = len(ind)
        state = 0
        for k in range(one_len):
            one_label = target_mask[ind[k, 0], ind[k, 1], ind[k, 2]]
            if one_label == 1:
                state += 1

        thres = <int> (threshold * one_len)
        if state > 0:
            if state >= thres:
                num_fiber += 1
                all_length += c_length(one_str, weighting)

    return num_fiber, all_length


def fiber_density_map(list streamlines, double [:, ::1] affine,
                      unsigned short [:] space):
    cdef int i, k, one_len
    cdef double [:,::1] one_str
    cdef unsigned short [:,::1]  ind
    cdef double [:,::1] lin_T
    cdef double [:] offset
    cdef unsigned short [:, :, :] density 

    density = np.zeros((space[0], space[1], space[2]), dtype=np.uint16)
    lin_T, offset = _mapping_to_voxel(affine)
    for i in range(len(streamlines)):
        one_str = streamlines[i]
        ind = _to_voxel_coordinates(one_str, lin_T, offset)
        one_len = len(ind)
        for k in range(one_len):
            density[ind[k, 0], ind[k, 1], ind[k, 2]] += 1

    return density
