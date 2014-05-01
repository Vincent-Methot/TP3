#!/usr/bin/python
# -*- coding: utf-8 -*-

"""TP3, question 4 (extra), IMN530. Par Jérémie Fouquet et Vincent Méthot"""


import numpy as np
import nibabel as nib
from dipy.core.gradients import gradient_table
from dipy.reconst.csdeconv import auto_response
from dipy.reconst.csdeconv import ConstrainedSphericalDeconvModel
from dipy.reconst.peaks import peaks_from_model
from dipy.data import get_sphere
from dipy.segment.mask import median_otsu

# Constrained Spherical Harmonics Reconstruction

dmri = nib.load('Data/dmri.nii')
data = dmri.get_data()
grad_dir = np.ndfromtxt('Data/gradient_directions_b-values.txt')
gtab = gradient_table(grad_dir[:, 3], grad_dir[:, :3])

response, ratio = auto_response(gtab, data, roi_radius=10, fa_thr=0.7)
csd_model = ConstrainedSphericalDeconvModel(gtab, response)
sphere = get_sphere('symmetric724')
csd_peaks = peaks_from_model(model=csd_model,
                             data=data,
                             sphere=sphere,
                             relative_peak_threshold=.5,
                             min_separation_angle=25,
                             parallel=True)

# Tracking with EuDX

from dipy.tracking.eudx import EuDX
eu = EuDX(csd_peaks.gfa,
          csd_peaks.peak_indices[..., 0],
          seeds=10000,
          odf_vertices=sphere.vertices,
          a_low=0.2)
csa_streamlines = [streamline for streamline in eu]
hdr = nib.trackvis.empty_header()
hdr['voxel_size'] = (2., 2., 2.)
hdr['voxel_order'] = 'LAS'
hdr['dim'] = csd_peaks.gfa.shape[:3]
csa_streamlines_trk = ((sl, None, None) for sl in csa_streamlines)
hdr['n_count'] = 10000
csa_sl_fname = 'Data/csd_streamline.trk'
nib.trackvis.write(csa_sl_fname, csa_streamlines_trk, hdr, points_space='voxel')
