import pickle as pkl
import numpy as np
from scipy import special as spec

import higrid.Microphone as mc


def shd_nm(channels, n, m):
    """
    Calculate the n-the order, m-th degree SHD coefficients for the given Eigenmike em32 channels

    :param channels: 32 channels of audio from Eigenmike em32
    :param n: SHD coefficient order
    :param m: SHD coefficient degree
    :return: Non-equalised SHD coefficients (p_nm) for the given em32 recording
    """
    em32 = mc.EigenmikeEM32()
    estr = em32.returnAsStruct()
    wts = estr['weights']
    ths = estr['thetas']
    phs = estr['phis']
    pnm = np.zeros(np.shape(channels[0])) * 1j
    for ind in range(32):
        cq = channels[ind]
        wq = wts[ind]
        tq = ths[ind]
        pq = phs[ind]
        Ynm = spec.sph_harm(m, n, pq, tq)
        pnm += wq * cq * np.conj(Ynm)
    return pnm


def shd_all(channels, Nmax = 4):
    """
    Calculate all SHD coefficients for the given em32 recording

    :param channels: 32 channels of audio from Eigenmike em32
    :param Nmax: Maximum SHD order to be calculated (default = 4)
    :return: List including (Nmax+1)^2 (complex-valued) SHD coefficients
    """
    Pnm = []
    for n in range(Nmax + 1):
        for m in range(-n, n + 1):
            pnm = shd_nm(channels, n, m)
            Pnm.append(pnm)
    return Pnm


def selectdeclevel(freq):
    """
    Calculate an appropriate SHD order for the given frequency

    :param freq: Frequency (in Hz)
    :return: Decomposition order

    Note: This is not used since we are looking at frequencies above 2607 Hz only. Hence Nmax=4 is used
    """
    if freq < 652.0:
        Ndec = 1
    elif freq < 1303.0:
        Ndec = 2
    elif freq < 2607.0:
        Ndec = 3
    else:
        Ndec = 4
    return Ndec


def getYnmtlookup(filepath):
    """
    Load the pre-computed spherical harmonics matrices

    :param filepath: Location of the pickled files
    :return: Pre-computed spherical harmonics matrices
    """
    f = open(filepath, 'rb')
    Ynmt = pkl.load(f)
    f.close()
    return Ynmt