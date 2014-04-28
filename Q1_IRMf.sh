#!/bin/bash

# Script de reconstruction des zones d'activation les plus importantes
# pour l'image Data/IRMf.nii. Utilise la série temporelle idéalle du
# fichier Data/ideal.txt

# Les opérations suivantes nécessitent les logiciels AFNI et FSL

##############################################################################

# Partie a - Étapes de reconstruction

# Extraction du cerveau (fsl-bet). Permet d'obenir une image plus petite
# et de minimiser le temps de calcul des opérations subséquentes.

bet Data/fmri.nii Data/fmri_bet.nii.gz -F -f 0.6

# Lissage spatial. Convolue l'image fmri.nii à chaque temps avec une
# gaussienne de FWHM = 6 mm.

3dmerge -1blur_fwhm 6 -doall -prefix Data/fmri_smooth.nii.gz \
	Data/fmri_bet.nii.gz

# Filtre passe-bas. La fréquence seuil est de 0.06 Hz. Puisque la tâche a une
# fréquence fondamentale de 1/50 sec = 0.02 Hz, cette opération devrait
# réduire le bruit sans toucher au signal cible.

3dBandpass -prefix Data/fmri_lowpass.nii.gz 0 0.06 Data/fmri_smooth.nii.gz

# Corrélation de chaque série temporelle avec le signal idéal
# (fichier ideal.txt). Donne une carte des coefficients de corrélation à
# chaqu'un des voxels

3dfim+ -input Data/fmri_lowpass.nii.gz -ideal_file Data/ideal.txt \
	-out 'Correlation' -bucket Data/fmri_corr.nii.gz

##############################################################################

# Partie b - Segmentation et étiquetage des zones d'activation

##############################################################################

# Partie c - Création de figures