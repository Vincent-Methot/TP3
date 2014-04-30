#! /bin/bash
# Q3_Fusion.sh

# Ce script reproduit les étapes importantes de recalage de la question 3 a. 
# Cependant, , ces étapes sont souvent longues et exigeantes en ressources 
# computaionnelles. Il n'est donc pas recommandé de lancer le script en une 
# seule fois, à moins de l'avoir bien acompris et analysé au préalable.


# 1. Preprocessing
# Extraction du cerveau dans le b0 (essentiel pour un bon recalage)
mrconvert -coord 3 0 Data/dmri.nii Data/b0.nii.gz
fsl5.0-bet Data/b0.nii.gz Data/b0_bet.nii.gz -R -f 0.25 -m
# Extraction du cerveau dans la t1 (essentiel pour un bon recalage)
fsl5.0-bet Data/t1.nii Data/t1_bet.nii.gz -R -m
# Extraction du cerveau dans toute la série fonctionnelle (essentiel)
bet Data/fmri.nii Data/fmri_bet.nii.gz -F -f 0.6
# Calcul de la moyenne temporelle de l'IRMf
fsl5.0-fslmaths Data/fmri_bet.nii.gz -Tmean Data/fmri_bet_tmean.nii.gz



# 2. Fusion de l'IRMd sur la t1
# Détermination de la transformation avec le b0_bet et le t1_bet
ANTS 3 -m CC[Data/t1_bet.nii.gz, Data/b0_bet.nii.gz, 1, 4] -o \
Data/b0_trans/b0_to_t1.nii.gz -t Elast[3] -i 30x20x10 -r Gauss[0.5,3]
# Application de la transformation sur tous les b-value de dmri
fsl5.0-fslsplit Data/dmri.nii.gz Data/dmri
for b0Im in Data/dmri00*
do
    WarpImageMultiTransform 3 b0Im b0Im -R Data/t1_bet.nii.gz \
    Data/b0_trans/b0_to_t1Warp.nii.gz Data/b0_tans/b0_to_t1Affine.txt
done
fsl5.0-fslmerge -t Data/dmri_to_t1.nii.gz Data/dmri00*
rm Data/dmri00*
# Ensuite, nous avons calculé les tenseurs pour l'IRMd recalée grâce à notre 
# fonction python Q2_IRMd.tenseur. Nous avons sauvegardé les tenseurs en format
# nifti, ce qui nous a permis de les ouvrir dans le Fibernavigator.



# 3. Fusion de l'IRMf sur la t1
# Détermination de la transformation avec le fmri_bet_tmean et le t1_bet
ANTS 3 -m CC[Data/t1_bet.nii.gz, Data/fmri_bet_tmean.nii.gz, 1, 4] \
-o Data/fmri_trans/fmri_to_t1.nii.gz -t Elast[3] -i 30x20x10 -r Gauss[0.5,3]
# Application de la transformation sur les blobs d'activation. Le --use-NN 
# assure que les étiquettes des blobs ne se verront pas attribuer des valeurs
# d'étiquette d'autres blobs (c'est un nearest neighbour).
WarpImageMultiTransform 3 Data/fmri_clusters.nii.gz \
Data/fmri_clusters_to_t1.nii.gz -R t1_bet.nii.gz \
Data/fmri_trans/fmri_to_t1Warp.nii.gz Data/fmri_trans/fmri_to_t1Affine.txt \
--use-NN

