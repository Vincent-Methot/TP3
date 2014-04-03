#!/usr/bin/python
# -*- coding: utf-8 -*-

"""TP3, question 2, IMN530"""

import nibabel as nib
import numpy as np
from dipy.core.ndindex import ndindex


def tenseur(dmri, gtab):

    dmri = nib.load(dmri)
    data = dmri.get_data()

    print 'Forme de dmri:', data.shape
    gtab = np.ndfromtxt(gtab)[1:]
    S0 = data[..., 0]
    S = data[..., 1:]

    B = np.array([ gtab[:,0]**2, gtab[:,0]*gtab[:,1], gtab[:,0]*gtab[:,2], gtab[:,1]**2, gtab[:,1]*gtab[:,2], gtab[:,2]**2 ])

    for index in ndindex(data.shape[:3]):
        X = -(1 / gtab[:, 3].astype(float)) * np.log( S[index].astype(float) / S0[index].astype(float) )



def compAdcAndFa(tensMat):
    """ Fonction qui calcule l'ADC et la fa d'une matrice 3D comprenant des
    tenseurs de diffusion à chaque index."""

    adcMap = np.zeros(tensMat.shape[:3])
    faMap = np.zeros(tensMat.shape[:3])

    for idx in ndindex(tensMat.shape[:3]):
        dLin = tensMat[idx]
        eigv = compLinDTensorEigval(dLin)
        adcMap[idx] = eigv.sum() / 3
        mEigv = eigv.mean()
        faMap[idx] = np.sqrt(3/2*((eigv-eigv.mean())**2).sum()/(eigv**2).sum())

    return adcMap, faMap



def compDiffTenEigval(dLin):
    dLin2MatIdx = np.array([[True, False, False],
                            [True, True, False],
                            [True, True, True]])
    d = np.zeros(3, 3)
    d[dLin2MatIdx] = dLin
    eigv = np.linalg.eigvalsh(d)

    return eigv