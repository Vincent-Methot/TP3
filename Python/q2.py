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
            X = -(1 / gtab[:, 3].astype(float)) * ( np.log( S[index] ) - np.log( S0[index] ) )
            tenseur[index] = np.dot(np.linalg.pinv(B), X)

    return tenseur


def compAdcAndFa(tensMat):
    """ Fonction qui calcule l'ADC et la fa d'une matrice 3D comprenant des
    tenseurs de diffusion à chaque index."""

    adcMap = np.zeros(tensMat.shape[:3])
    faMap = np.zeros(tensMat.shape[:3])

    for idx in ndindex(tensMat.shape[:3]):
        dLin = tensMat[idx]
        eigv = compLinDTensorEigval(dLin)
        adcMap[idx] = eigv.sum() / 3
        faMap[idx] = np.sqrt(3/2*((eigv-eigv.mean())**2).sum()/(eigv**2).sum())

    return adcMap, faMap


def compLinDTensorEigval(dLin):
    pdb.set_trace()
    dLin2MatIdx = np.array([[True, True, True],
                            [False, True, True],
                            [False, False, True]])
    dt = np.zeros([3, 3])
    dt[dLin2MatIdx] = dLin
    eigv = np.linalg.eigvalsh(dt, UPLO='U')

    return eigv