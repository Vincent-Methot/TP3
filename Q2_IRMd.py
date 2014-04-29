#!/usr/bin/python
# -*- coding: utf-8 -*-

"""TP3, question 2, IMN530. Par Jérémie Fouquet et Vincent Méthot"""

import nibabel as nib
import numpy as np
import scipy as sp
from dipy.core.ndindex import ndindex
from dipy.segment.mask import median_otsu
import pdb


def tenseur(dmri, gtab, bMaskSource=None):
    """Estimation des tenseurs par méthode des moindres carrés.

    Paramètres
    ----------
    dmri: Fichier nifti contenant une image d'IRM de diffusion
    gtab: Table des gradients utilisée
    bMaskSource: Fichier nifti contenant un masque binaire du cerveau.

    Retour
    ------
    tenseur: nparray contenant les tenseurs à chaque voxel

    Example
    -------
    >> tenseur = q2.tenseur('Data/dmri.nii',
                 '../Data/gradient_directions_b-values.txt')"""

    # Chargement des données et déclaration de variables
    dmri = nib.load(dmri)
    data = dmri.get_data()

    gtab = np.ndfromtxt(gtab)[1:]
    S0 = data[..., 0]
    S = data[..., 1:]

    B = np.array([gtab[:, 0] ** 2, gtab[:, 0] * gtab[:, 1],
                  gtab[:, 0] * gtab[:, 2], gtab[:, 1] ** 2,
                  gtab[:, 1] * gtab[:, 2], gtab[:, 2] ** 2]).T
    tenseur = np.empty(data.shape[:3] + (6,))

    # Définition des indices d'itération en fonction du masque de cerveau
    if (bMaskSource is None) | (bMaskSource is False):
        itIdx = ndindex(data.shape[:3])
    else:
        bMask = openimage(bMaskSource)
        itIdx = bMask.nonzero()
        itIdx = list(zip(itIdx[0], itIdx[1], itIdx[2]))

    # Calcul des tenseurs en tous les indices définis précédemment
    for index in itIdx:
        if S0[index] == 0:
            tenseur[index] = np.zeros(6)
        else:
            X = -((1 / gtab[:, 3].astype(float)) *
                (np.log(S[index].astype(float) / S0[index].astype(float))))
            tenseur[index] = np.dot(np.linalg.pinv(B), X)

    tenseur[np.isinf(tenseur) | np.isnan(tenseur)] = 0
    return tenseur.astype('f32')


def compAdcAndFa(tensMat, bMaskSource=None):
    """ Calcul de l'ADC et de la FA d'une matrice 3D comprenant des
    tenseurs de diffusion à chaque index.

    Paramètres
    ----------
    tensMat:nparray. Matrice MxNxPx6 comprenant les tenseurs
    bMaskSource: Fichier nifti contenant un masque pour lequel l'ADC et la FA
    seront calculées (généralement un masque de cerveau).

    Retour
    ------
    adcMap: nparray. Matrice MxNxP de l'ADC a chaque voxel
    faMap: nparray. Matrice MxNxP de la FA a chaque voxel"""

    # Initialisation des variables de résultat
    adcMap = np.zeros(tensMat.shape[:3])
    faMap = np.zeros(tensMat.shape[:3])

    # Définition des indices d'itération en fonction du masque de cerveau
    if (bMaskSource is None) | (bMaskSource is False):
        itIdx = ndindex(tensMat.shape[:3])
    else:
        bMask = openimage(bMaskSource)
        itIdx = bMask.nonzero()
        itIdx = list(zip(itIdx[0], itIdx[1], itIdx[2]))

    # Calcul des cartes de FA et d'ADC
    for idx in itIdx:
        dLin = tensMat[idx]
        eigv = compLinDTensorEigv(dLin)
        adcMap[idx] = eigv.sum() / 3
        faMap[idx] = np.sqrt(3 / 2 * ((eigv - eigv.mean()) ** 2).sum() /
                     (eigv ** 2).sum())

    return adcMap.astype('f32'), faMap.astype('f32')


def compLinDTensorEigv(dLin, compEigVec=False):
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

    if compEigVec:
        return np.linalg.eigh(dt, UPLO='U')
    else:
        return np.linalg.eigvalsh(dt, UPLO='U')


def tracking(tensMat, trackStep=0.5, nSeed=10000, faTh=0.15, maxAngle=np.pi/3,
             splineOrder=1, wmMaskSource=None, bMaskSource=None, fa=None,
             saveTracksFname=None, trackHdr=None, verbose=True,
             dataToInterp='tens'):
    """Tracking déterministe de fibre dans la matrice de tenseurs tensMat, dans
    un masque de matière blanche wmMaskSource (fichier nifti).

    Paramètres
    ----------
    tensMat: nparray MxNxPx6 contenant les tenseurs de diffusion
    trackStep: float. Facteur multiplicatif qui dicte la longueur des pas.
    nSeed: int. Nombre de seeds placées aléatoirement dans la matière blanche
    faTh: float. Threshold sur la FA pour déterminer où faire le tracking.
    maxAngl: float. Angle entre deux directions de tracking successives
        au-dessus duquel le tracking est arrêté.
    wmMaskSource: Fichier nifti contenant un masque binaire de la matière
        blanche. Si None, le masque est calculé à partir des param. restants.
    bMaskSource: Fichier nifti contenant un masque binaire MxNxP du cerveau
        (utilisé seulement si wmMaskSource=None).
    fa: nparray MxNxP contenant la FA (utilisé seulement si wmMaskSource=None)
    saveTracksFname: String. Nom d'un fichier .trk dans lequel les tracks
        trouvées vont être sauvegardées
    trackHdr: nifti header (nibabel). Header duquel des infos sont prises pour
    sauvegarder les tracks. Si None, pour l'instant hdr['voxel_size']=(1,1,1).

    Retour
    ------
    allPts: Les points des streamlines trouvés pour toutes les seeds, sous la
    forme d'une liste de longueur nSeed contenant des nparrray nbOfPts x 3.

    Exemple
    -------
    >> allPts = Q2_IRMd.tracking(tenseur, bMaskSource='Data/b0_bet_mask.nii.gz',
                                 fa='Data/FA.nii.gz', verbose=False,
                                 saveTracksFname='Data/tracks.trk')
    """

    # 1. Détermination du masque de la matiere blanche
    if wmMaskSource:
        wmMask = openimage(wmMaskSource)
    else:
        if not(fa):
            dum, fa = compAdcAndFa(tensMat, bMaskSource)
            del dum
            # Il resterait à implémenter une façon de calculer un masque du
            # cerveau si bMaskSource == None.
        wmMask = segmentwhitematter(bMaskSource, fa, faTh)

    # 2. Calcul des informations nécessaires en fonction de ce que l'on veut
    # interpoler. À voir si le temps le permet!

    # 3. Détermination aléatoire des coordonnées des seeds
    seedPts = wmMask.nonzero()
    np.random.seed(1) # Seulement utilisé pour des fins d'optimisation. À enlever par la suite.
    seedIdx = np.random.permutation(seedPts[0].size)
    seedIdx = list(seedIdx[list(range(nSeed))])
    seedPts = list(zip(seedPts[0][seedIdx], seedPts[1][seedIdx],
                   seedPts[2][seedIdx]))

    # 4. Tracking
    nextTens = np.zeros(6,)
    minAngleCos = np.cos(maxAngle)
    allPts = []
    for iSeed in range(nSeed):
        # Affichage du no de seed
        if verbose:
            print(("{}{}".format("Seed ", iSeed)))
        # Initialisation de la direction de tracking
        actualPt = seedPts[iSeed]
        ptList = [np.array(actualPt)]
        actualEva, actualEve = compLinDTensorEigv(tensMat[actualPt[0],
                                actualPt[1], actualPt[2], :], compEigVec=True)
        maxEvaIdx = actualEva.argmax()
        # Boucle sur les 2 directions possibles à partir de la seed
        for iDir in range(2):
            actualDir = actualEve[:, maxEvaIdx]
            if iDir == 1:
                actualDir = -actualDir
            continueTrack = True
            outOfMask = False
            maxAngleReached = False
            while continueTrack:
                # On avance d'un pas
                nextPt = actualPt + actualDir * trackStep
                # On regarde si on sort du masque
                if not(wmMask[tuple(np.round(nextPt).astype(int))]):
                    outOfMask = True
                    break
                # On évalue la nouvelle direction
                if splineOrder == 'nearest':
                    interpPt = np.round(nextPt).astype(int)
                    nextTens = tensMat[interpPt[0], interpPt[1],
                                       interpPt[2], ...]
                else:
                    interpPt = np.array([nextPt])
                    for iTens in range(6):
                        nextTens[iTens] = sp.ndimage.map_coordinates(
                                            tensMat[..., iTens], interpPt.T,
                                            order=splineOrder)
                nextEva, nextEve = compLinDTensorEigv(tensMat[nextPt[0],
                                    nextPt[1], nextPt[2], :],
                                                          compEigVec=True)
                maxEvaIdx = nextEva.argmax()
                nextDir = nextEve[:, maxEvaIdx]
                # On vérifie l'angle entre l'ancienne et la nouvelle direction
                if np.dot(nextDir, actualDir) < 0:
                    nextDir = -nextDir
                if np.dot(nextDir, actualDir) < minAngleCos:
                    maxAngleReached = True
                    break
                continueTrack = not(outOfMask | maxAngleReached)
                # On met à jour et sauvegarde le point de tracking
                actualPt = nextPt
                actualDir = nextDir
                if iDir == 0:
                    ptList = ptList + [np.array(actualPt)]
                else:
                    ptList = [np.array(actualPt)] + ptList

        # Stockage des points trouvés pour la seed
        allPts = allPts + [np.array(ptList)]

    # Save trakcs in .trk file if file name is provided (par Maxime Descoteaux)
    if saveTracksFname:
        streamlines_trk = ((sl, None, None) for sl in allPts)
        # Construct header
        print('Saving fibers in trk format...')
        hdr = nib.trackvis.empty_header()
        if trackHdr:
            hdr['voxel_size'] = trackHdr.get_zooms()[:3]
        else:
            hdr['voxel_size'] = (1, 1, 1)
        hdr['voxel_order'] = 'LAS'
        hdr['dim'] = tensMat.shape[:3]
        hdr['n_count'] = len(allPts)
        nib.trackvis.write(saveTracksFname, streamlines_trk, hdr,
                            points_space='voxel')

    return allPts


def segmentwhitematter(bMaskSource, fa, faTh=0.15):
    """ Détermine le masque de la matière blanche à partir d'un masque
    du cerveau (fichier nifti) et d'un seuil sur la FA. """

    # Ouverture du fichier du masque du cerveau
    bMask = openimage(bMaskSource)

    # Seuil sur la fa
    #faTh = np.median(fa[bMask.astype(bool)])
    faMask = (fa > faTh)

    # Intersection des deux masques
    wmMask = bMask & faMask

    return wmMask


def compmainevec(tensMat, bMaskSource=None):
    """Calcul du vecteur propre avec la plus grande valeur propre à chaque
    position du champ de tenseur tensMat"""

    # Initialisation des variables de résultat
    mainVec = np.zeros(tensMat.shape[:3] + (3,))
    allEve = np.zeros(tensMat.shape[:3] + (3,3))
    allEva = np.zeros(tensMat.shape[:3] + (3,))

    # Définition des indices d'itération en fonction du masque de cerveau
    if (bMaskSource is None) | (bMaskSource is False):
        itIdx = ndindex(tensMat.shape[:3])
    else:
        bMask = openimage(bMaskSource)
        itIdx = bMask.nonzero()
        itIdx = list(zip(itIdx[0], itIdx[1], itIdx[2]))

    # Calcul du vecteur propre principal
    for idx in itIdx:
        dLin = tensMat[idx]
        eva, eve = compLinDTensorEigv(dLin, compEigVec=True)
        mainEvIdx = eva.argmax()
        mainVec[idx] = eve[:, mainEvIdx]
        allEve = eve
        allEva = eva

    return mainVec, allEve, allEva


def openimage(im):
    """Fonction permettant de transformer mask en nparray. Mask peut être un
    nparray, un fichier nifti, ou une liste."""

    if isinstance(im, str):
        data = (nib.load(im)).get_data()
    else:
        data = np.array(im)

    return data
