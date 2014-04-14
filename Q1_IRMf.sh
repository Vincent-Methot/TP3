#!/bin/bash

# Script de reconstruction des zones d'activation les plus importantes
# pour l'image Data/IRMf. Utilise les 'ideal time course' du fichier
# Data/ideal.txt

# Extraction du cerveau (fsl-bet)
bet fmri.nii fmri_bet.nii.gz -F

# Lissage spatial
3dmerge -1blur_fwhm 6 -doall -prefix fmri_smooth.nii.gz fmri_bet.nii.gz

# Filtre passe-bas
3dBandpass -prefix fmri_lowpass.nii.gz 0 .06 fmri_smooth.nii.gz
