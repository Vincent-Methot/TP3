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
