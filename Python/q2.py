#!/usr/bin/python
# -*- coding: utf-8 -*-

"""TP3, question 2, IMN530"""

import nibabel as nib
import numpy as np
from dipy.core.ndindex import ndindex
import pdb


def tenseur(dmri, gtab):
    """Estimation des tenseurs par méthode des moindres carrés.

    Paramètres
    ----------
    dmri: Fichier nifti contenant une image d'IRM de diffusion
    gtab: Table des gradients utilisée

    Retour
    ------
    tenseur: nparray contenant les tenseurs à chaque voxel"""

    dmri = nib.load(dmri)
    data = dmri.get_data()

    print 'Forme de dmri:', data.shape
    gtab = np.ndfromtxt(gtab)[1:]
    S0 = data[..., 0]
    S = data[..., 1:]

    B = np.array([ gtab[:,0]**2, gtab[:,0]*gtab[:,1], gtab[:,0]*gtab[:,2], gtab[:,1]**2, gtab[:,1]*gtab[:,2], gtab[:,2]**2 ]).T
    tenseur = np.empty(data.shape[:3] + (6,))

    for index in ndindex(data.shape[:3]):
        if S0[index] == 0:
            tenseur[index] = np.zeros(6)
        else:
            X = -(1 / gtab[:, 3].astype(float)) * ( np.log( S[index].astype(float) / S0[index].astype(float) ) )
            tenseur[index] = np.dot( np.linalg.pinv(B), X )

    tenseur[np.isinf(tenseur) | np.isnan(tenseur)] = 0
    return tenseur


def compAdcAndFa(tensMat):
    """ Calcul de l'ADC et de la FA d'une matrice 3D comprenant des
    tenseurs de diffusion à chaque index.

    Paramètres
    ----------
    tensMat:nparray. Matrice MxNxPx6 comprenant les tenseurs

    Retour
    ------
    adcMap: nparray. Matrice MxNxP de l'ADC a chaque voxel
    faMap: nparray. Matrice MxNxP de la FA a chaque voxel"""


    adcMap = np.zeros(tensMat.shape[:3])
    faMap = np.zeros(tensMat.shape[:3])

    for idx in ndindex(tensMat.shape[:3]):
        dLin = tensMat[idx]
        eigv = compLinDTensorEigval(dLin)
        adcMap[idx] = eigv.sum() / 3
        faMap[idx] = np.sqrt(3/2*((eigv-eigv.mean())**2).sum()/(eigv**2).sum())

    return adcMap, faMap


def compLinDTensorEigval(dLin, compEigVec=False):
    """Calcul des valeurs propres d'un tenseur symétrique T 3x3 exprimé sous
    la forme d'un vecteur V 6x1 tel que V[0]=T[0,0], V[1]=T[0,1], V[2]=T[0,2],
    V[3]=T[1,1], V[4]=T[1,2], V[5]=T[2,2]

    Paramètres
    ----------
    dLin:nparray. Tenseur 1x6

    Retour
    ------
    eigv: nparray. Valeurs propres du tenseurs
    """
    dLin2MatIdx = np.array([[True, True, True],
                            [False, True, True],
                            [False, False, True]])
    dt = np.zeros([3, 3])
    dt[dLin2MatIdx] = dLin

    if ~compEigVec:
        return np.linalg.eigvalsh(dt, UPLO='U')
    else:
        return np.linalg.eigh(dt, UPLO='U')


def tracking(tensMat):
    """Tracking déterministe de fibre dans la matrice de tenseurs tensMat."""

    # Détermination du masque de la matière blanche
