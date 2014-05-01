#!/bin/bash

# Script de reconstruction des zones d'activation les plus importantes
# pour l'image Data/IRMf.nii. Utilise la série temporelle idéalle du
# fichier Data/ideal.txt

# Les opérations suivantes nécessitent les logiciels AFNI et FSL

##############################################################################

# Partie a - Étapes de reconstruction

# Extraction du cerveau (fsl-bet). Permet d'obenir une image plus petite
# et de minimiser le temps de calcul des opérations subséquentes.

fsl5.0-bet Data/fmri.nii Data/fmri_bet.nii.gz -F -f 0.6

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
3dclust -savemask Data/fmri_clusters.nii.gz -1Dformat -nosum -1dindex 0 \
-1tindex 0 -1noneg -2thresh -0.5149 0.5149 -dxyz=1 1.01 50 Data/fmri_corr.nii.gz

# On obtient les différentes régions d'activation, portant les étiquettes:
#
# 1 - Aire visuelle droite
# 2 - Aire visuelle gauche
# 3 - Hippocampe gauche
# 4 - Hippocampe droit
# 5 - Cortex moteur droit
# 6 - Artéfact hors du cerveau
# 7 - Cortex moteur gauche

##############################################################################

# Partie c - Création de figures

# Calcul de l'image moyenne de fmri.nii pour des fins de visualisation dans
# le fibernavigator.

3dTstat -mean -prefix Data/fmri_mean.nii.gz Data/fmri_bet.nii.gz 