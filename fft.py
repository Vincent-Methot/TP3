#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Fait la fft de l'image fmri.nii (qui doit être dans le présent dossier
et retourne fmri_fft.nii"""

import nibabel as nib
import numpy as np

fmri = nib.load('fmri.nii')
data = fmri.get_data()
affine = fmri.get_affine()
hdr = fmri.get_header()
tr = hdr.get_zooms()[-1]
zooms = hdr.get_zooms()

freq = np.fft.rfftfreq(data.shape[-1], d=tr)
fourier = abs(np.fft.rfft(data))
df = freq[1] - freq[0]
zooms = (zooms[0], zooms[1], zooms[2], df)
hdr.set_zooms(zooms)

Fourier = nib.Nifti1Image(fourier, affine, hdr)
nib.save(Fourier, 'fmri_fft.nii')
