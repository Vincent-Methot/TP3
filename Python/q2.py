#!/usr/bin/python
# -*- coding: utf-8 -*-

"""TP3, question 2, IMN530"""

import nibabel as nib
import numpy as np
from dipy.core.ndindex import ndindex


def tenseur(dmri, gtab):
	dmri = nib.load('dmri.nii')
	data = dmri.get_data()

	print 'Forme de dmri.nii:', data.shape
	gtab = np.ndfromtxt('gradient_directions_b-values.txt')[1:]
	S0 = data[..., 0]
	S = data[..., 1:]

	for index in ndindex(data.shape[:3]):
		X = -(1 / gtab[:, 3]) * np.log( S[index] / S0[index] )
