import nibabel as nib
import numpy as np

dmri = nib.load('Data/dmri.nii')
data = dmri.get_data()
hdr = dmri.get_header()
affine = dmri.get_affine()
b0 = data[...,0]
B0 = nib.Nifti1Image(b0, affine())
nib.save(B0, 'b0.nii.gz')

# bet b0.nii.gz b0_bet.nii.gz -m

mask = nib.load('Data/b0_bet_mask.nii.gz').get_data()
dmri_brain = np.empty(data.shape)

# Pas la meilleure façon de faire, mais la plus simple...
for i in range(data.shape[-1]):
	dmri_brain[..., i] = data[..., i] * mask


