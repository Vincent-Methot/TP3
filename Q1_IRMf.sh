#!/bin/bash

# Script de reconstruction des zones d'activation les plus importantes
# pour l'image Data/IRMf. Utilise les 'ideal time course' du fichier
# Data/ideal.txt

# Nécessite les logiciels AFNI et FSL

# Extraction du cerveau (fsl-bet)
bet fmri.nii fmri_bet.nii.gz -F -f 0.6

# Lissage spatial
3dmerge -1blur_fwhm 6 -doall -prefix fmri_smooth.nii.gz fmri_bet.nii.gz

# Filtre passe-bas
3dBandpass -prefix fmri_lowpass.nii.gz 0 .06 fmri_smooth.nii.gz

# Corrélation
3dfim+ -input fmri_lowpass.nii.gz -ideal_file ideal.1D -out 'Correlation' -bucket fmri_corr.nii.gz
