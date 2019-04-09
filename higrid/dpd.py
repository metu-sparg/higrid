from collections import defaultdict
import numpy as np
from scipy import special as spec
from higrid.utils import sph_jnyn


def getWY(micstruct, Ndec):
    """
    Return the array (e.g. em32) specific cubature and the SHD matrices, W and Y.H

    :param micstruct: Dictionary containing microphone properties
    :param Ndec: SHD order to be used in the cubature
    :return: Array specific cubature and the SHD matrices, W and Y.H
    """
    thes = micstruct['thetas']
    phis = micstruct['phis']
    w = micstruct['weights']
    W = np.matrix(np.diag(w), dtype=float)
    Y = np.matrix(np.zeros((len(thes), (Ndec + 1) ** 2), dtype=complex))
    for ind in range(len(thes)):
        y = []
        for n in range(Ndec + 1):
            for m in range(-n, n + 1):
                Ynm = spec.sph_harm(m, n, phis[ind], thes[ind])
                y.append(Ynm)
        Y[ind, :] = y
    W = W / np.diag(Y * Y.H) / 2.
    return W, Y.H


def getBmat(micstruct, findmin, findmax, NFFT, Fs, Ndec):
    """
    Return the array (e.g. em32) specific response equalisation matrix, B

    :param micstruct: Dictionary containing microphone properties
    :param findmin: Index of minimum frequency (int)
    :param findmax: Index of maximum frequency (int)
    :param NFFT: FFT size
    :param Fs: Sampling rate (Hz)
    :param Ndec: SHD order
    :return: Response equalisation matrices, B, for each frequency index
    """
    ra = micstruct['radius']
    Bmat = defaultdict()
    for find in range(findmin, findmax):
        freq = float(find) * Fs / NFFT
        kra = 2 * np.pi * freq / 344.0 * ra
        jn, jnp, yn, ynp = sph_jnyn(Ndec, kra)
        # jn, jnp, yn, ynp = spec.sph_jnyn(Ndec, kra) # scipy 0.19.1
        hn = jn - 1j * yn
        hnp = jnp - 1j * ynp
        bnkra = jn - (jnp / hnp) * hn
        bval = []
        for ind in range(Ndec + 1):
            for jnd in range(2 * ind + 1):
                bval.append(bnkra[ind] * 4 * np.pi * (1j) ** ind)
        Bmat[find] = np.linalg.inv(np.matrix(np.diag(bval)))
    return Bmat

def getpvec(P, tind, find):
    """
    Return a single M-channel (e.g. 32 channel for em32) time-frequency bin

    :param P: List of STFTs of each channel
    :param tind: Time index
    :param find: Frequency index
    :return: Selected time frequency bin containing N (e.g. 32) channels
    """
    pvec = []
    for ind in range(len(P)):
        pvec.append(P[ind][tind, find])
    pvec = np.matrix(pvec)
    return pvec.T


def getanmval(pvec, B, Y, W):
    """
    Return the SHD for a single time-frequency bin

    :param pvec: Vector containing M-channel STFTs of a single time-frequency bin
    :param B: Response equalisation matrix
    :param Y: SHD matrix
    :param W: Cubature matrix
    :return: SHD for a single time-frequency bin; (N+1)^2 by 1
    """
    anm = B * Y * W * pvec
    return anm


def getAnm(P, mstr, Bmat, findmin, findmax, Ndec):
    """
    Return the (N+1)^2-element list containing SHDs of STFTs

    :param P: STFTs of the M channels of recordings
    :param mstr: Dict containing the microphone array properties
    :param Bmat: List of response equalisation matrices
    :param findmin: Index of minimum frequency (int)
    :param findmax: Index of maximum frequency (int)
    :param Ndec: SHD order
    :return: List of numpy matrices containing the SHDs of STFTs of array channels
    """
    A = []
    W, Y = getWY(mstr, Ndec)
    for ind in range((Ndec + 1) ** 2):
        A.append(np.zeros((P[0].shape[0], P[0].shape[1]), dtype=complex))
    for find in range(findmin, findmax):
        B = Bmat[find]
        for tind in range(P[0].shape[0]):
            pv = getpvec(P, tind, find)
            anm = getanmval(pv, B, Y, W)
            for snd in range((Ndec + 1) ** 2):
                A[snd][tind, find] = anm[snd]
    return A


def dpd(Anm, Ndec, find, tind, Jtau, Jnu, thr):
    """
    Direct Path Dominance (DPD) test

    :param Anm: (N+1)^2-element list containing SHDs of STFTs
    :param Ndec: SHD order
    :param find: Frequency index
    :param tind: Time index
    :param Jtau: Time averaging window size
    :param Jnu: Frequency averaging window size
    :param thr: DPD threshold
    :return: flag (1 if erank=1, 0 otherwise)

    Note:
    See the following paper for details of DPD
    Nadiri, O., and Rafaely, B. (2014). Localization of multiple speakers under high reverberation using a spherical
    microphone array and the direct-path dominance test. IEEE/ACM Trans. on Audio, Speech, and Lang. Process., 22(10),
    1494-1505.
    """
    anm = np.matrix(np.zeros(((Ndec + 1) ** 2, 1), dtype=complex))
    Ra = np.matrix(np.zeros(((Ndec + 1) ** 2, (Ndec + 1) ** 2), dtype=complex))

    for fi in range(find, find + Jnu):
        for ti in range(tind, tind + Jtau):
            for ind in range((Ndec + 1) ** 2):
                anm[ind, 0] = Anm[ind][ti, fi]
            Ra += (anm * anm.H)
    Ra = Ra / (Jtau * Jnu)
    S = np.linalg.svd(Ra, compute_uv=False)
    ratio = S[0] / S[1]
    if ratio > thr:
        flag = 1
    else:
        flag = 0
    return flag