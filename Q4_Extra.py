#!/usr/bin/python
# -*- coding: utf-8 -*-

"""TP3, question 4 (extra), IMN530. Par Jérémie Fouquet et Vincent Méthot"""

import numpy as np
import nibabel as nib
from dipy.data import fetch_stanford_hardi, read_stanford_hardi, get_sphere
from dipy.reconst.shm import CsaOdfModel, normalize_data
from dipy.reconst.peaks import peaks_from_model