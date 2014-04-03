#!/usr/bin/python
# -*- coding: utf-8 -*-

"""TP3, question 2, IMN530"""

import nibabel as nib
import numpy as np
from dipy.core.ndindex import ndindex


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