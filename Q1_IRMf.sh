#!/bin/bash

# Script de reconstruction des zones d'activation les plus importantes
# pour l'image Data/IRMf.nii. Utilise la série temporelle idéalle du
# fichier Data/ideal.txt

# Nécessite les logiciels AFNI et FSL

# Extraction du cerveau (fsl-bet)
bet Data/fmri.nii Data/fmri_bet.nii.gz -F -f 0.6

# Lissage spatial
3dmerge -1blur_fwhm 6 -doall -prefix Data/fmri_smooth.nii.gz Data/fmri_bet.nii.gz

# Filtre passe-bas
3dBandpass -prefix Data/fmri_lowpass.nii.gz 0 .06 Data/fmri_smooth.nii.gz

# Corrélation
3dfim+ -input Data/fmri_lowpass.nii.gz -ideal_file Data/ideal.txt -out 'Correlation' -bucket Data/fmri_corr.nii.gz
