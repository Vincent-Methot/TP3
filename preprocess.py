import nibabel as nib
import numpy as np
import Q2_IRMd

dmri = nib.load('Data/dmri.nii')
data = dmri.get_data()
hdr = dmri.get_header()
affine = dmri.get_affine()
data_dtype = dmri.get_data_dtype()
b0 = data[...,0]
B0 = nib.Nifti1Image(b0, affine)
nib.save(B0, 'Data/b0.nii.gz')

# Extraction du cerveau dans fsl
# bet Data/b0.nii.gz Data/b0_bet.nii.gz -m

mask = nib.load('Data/b0_bet_mask.nii.gz').get_data()
data_brain = np.empty(data.shape)

# Pas la meilleure façon de faire, mais la plus simple...
for i in range(data.shape[-1]):
	data_brain[..., i] = data[..., i] * mask

dmri_brain = nib.Nifti1Image(data_brain.astype(data_dtype), affine)
nib.save(dmri_brain.set_data_dtype(data_dtype), 'Data/dmri_brain.nii.gz')

# Création de dmri_petit et dmri_mini

data_petit = data[30:95, 20:100, 0:55, :]
dmri_petit = nib.Nifti1Image(data_petit, affine)
nib.save(dmri_petit, 'Data/dmri_petit.nii.gz')
data_mini = data[50:70, 50:70, 20:40, :]
dmri_mini = nib.Nifti1Image(data_mini, affine)
nib.save(dmri_mini, 'Data/dmri_mini.nii.gz')

# Calcul des tenseurs, FA et ADC

tenseur_mini = Q2_IRMd.tenseur('Data/dmri_mini.nii.gz','Data/gradient_directions_b-values.txt')
tenseur_petit = Q2_IRMd.tenseur('Data/dmri_petit.nii.gz', 'Data/gradient_directions_b-values.txt')
tenseur = Q2_IRMd.tenseur('Data/dmri.nii', 'Data/gradient_directions_b-values.txt', 'Data/b0_bet_mask.nii.gz')

nib.save(nib.Nifti1Image(tenseur_mini, affine), 'Data/tenseur_mini.nii.gz')
nib.save(nib.Nifti1Image(tenseur_petit, affine), 'Data/tenseur_petit.nii.gz')
nib.save(nib.Nifti1Image(tenseur, affine), 'Data/tenseur.nii.gz')

ADC_mini, FA_mini = Q2_IRMd.compAdcAndFa(tenseur_mini)
ADC_petit, FA_petit = Q2_IRMd.compAdcAndFa(tenseur_petit)
ADC, FA = Q2_IRMd.compAdcAndFa(tenseur, 'Data/b0_bet_mask.nii.gz')

nib.save(nib.Nifti1Image(FA_mini, affine), 'Data/FA_mini.nii.gz')
nib.save(nib.Nifti1Image(FA_petit, affine), 'Data/FA_petit.nii.gz')
nib.save(nib.Nifti1Image(FA, affine), 'Data/FA.nii.gz')
nib.save(nib.Nifti1Image(ADC_mini, affine), 'Data/ADC_mini.nii.gz')
nib.save(nib.Nifti1Image(ADC_petit, affine), 'Data/ADC_petit.nii.gz')
nib.save(nib.Nifti1Image(ADC, affine), 'Data/ADC.nii.gz')

# Tracking...

allPts = Q2_IRMd.tracking(tenseur, bMaskSource='Data/b0_bet_mask.nii.gz',
						  fa='Data/FA.nii.gz', verbose=False,
						  saveTracksFname='Data/tracks.trk')