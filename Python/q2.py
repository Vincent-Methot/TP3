#!/usr/bin/python
# -*- coding: utf-8 -*-

"""TP3, question 3, IMN530"""

import nibabel as nib
import numpy as np

dmri = nib.load('dmri.nii')
data = dmri.get_data()
hdr = dmri.get_header()

print 'Forme de dmri.nii:', data.shape

fichier = open('gradient_directions_b-values.txt')
gradients = fichier.read()
gradients = gradients.split('\n')[:-1]
gradients = [s.split('\t') for s in gradients]
gradients = np.array(gradients)[:,:-1]