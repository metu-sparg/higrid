import numpy as np
import healpy as hp
from collections import defaultdict
import madmom as mm
from peakutils import indexes, baseline, scale
from random import shuffle
from scipy.ndimage.filters import median_filter, gaussian_filter
import tqdm

from higrid.dpd import getBmat, getAnm, dpd
from higrid.shd import shd_all
from higrid.treeutils import parent, children
from higrid.utils import node2vec, selectsome, cart2sph, selectbinindx, histtotree, sph_jnyn
from higrid.Microphone import EigenmikeEM32

import higrid.TreeClusteri as tc


def srpd(mvec, k, ra, Nmax, w, V):

    """
    Calculate the Steered Response Power Density (SRPD)

    :param mvec: SHD coefficients for the TF bin to be analysed
    :param k: Wave number (2*pi*f/c)
    :param ra: Radius of the microphone array
    :param Nmax: Maximum SHD order to be used
    :param w: Diagonal eigenvalue matrix
    :param V: Reduced eigenvector matrix
    :return: SRPD for the given pixel
    """
    assert np.size(mvec) == (Nmax + 1) ** 2
    V = V[0:(Nmax + 1) ** 2, 0:(Nmax + 1) ** 2]
    w = w[0:(Nmax + 1) ** 2]
    kra = k * ra
    jn, jnp, yn, ynp = sph_jnyn(Nmax, kra)
    # jn, jnp, yn, ynp = spec.sph_jnyn(Nmax, kra)
    hn = jn - 1j * yn
    hnp = jnp - 1j * ynp
    bnkra = jn - (jnp / hnp) * hn
    b = []
    for n in range(Nmax + 1):
        for count in range(-n, n + 1):
            b.append(1 / (4 * np.pi * (1j) ** n * bnkra[n]))
    b = np.array(b)
    p = b * mvec
    B0 = np.conj(np.matrix(np.conj(p)) * V).T
    B0s = np.diag(w) * np.multiply(B0, np.conj(B0))
    srpval = B0s.sum()
    return srpval

def treeroott(mvec, k, ra, Nmax, pixeldict, rlevel=1):
    """
    Calculate  SRPD for the all of the 3x2^(2nside) pixels at the lowest resolution level

    :param mvec: SHD coefficients of the TF bin to be analysed
    :param k: Wave number (2*pi*f/c)
    :param ra: Radius of the microphone array
    :param Nmax: Maximum order of the SHD to be used
    :param pixeldict: Precalculated eigenvalue and eigenvector matrices (see loadpixbasis in utils.py)
    :param rlevel: Resolution level (default=1; corresponds to the lowest level with 12 pixels)
    :return: Dictionary with keys in the form (level, index) and values containing SRPD value of that pixel
    """
    tree = defaultdict()
    nleaf = hp.nside2npix(rlevel)
    for ind in range(nleaf):
        w, V = pixeldict[(rlevel-1, ind)]
        val = srpd(mvec, k, ra, Nmax, w, V)
        tree[(0, ind)] = val
    return tree

def treeEntropy(tr):
    """
    Return spatial entropy of the representation

    :param tr: Tree representation of the SRPD map (as a defaultdict object)
    :return: Total spatial entropy of the representation

    Note: We are using the Cython version of this function. Please see cTreei.pyx
    """
    ln = len(tr)
    pw = np.zeros(ln, dtype=float)
    areas = np.zeros(ln, dtype=float)
    ind = 0
    for nd in tr:
        val = tr[nd]
        pw[ind] = np.abs(val)
        areas[ind] = ((2 ** (2 * nd[0]) * 12)) / (4 * np.pi)
        ind += 1
    pi = pw / np.sum(pw)
    etrpy = -np.sum(pi * np.log(pi * areas))  # / np.log(ln)
    return etrpy



def selecttreet(tr, er, level, idx, mvec, k, ra, Nmax, pixeldict):
    """
    Expand the multiresolution HiGRID map representation based on its total spatial entropy

    :param tr: Tree representation of the SRPD map (as a defaultdict object)
    :param er: Total spatial entropy of the representation before increasing resolution
    :param level: Resolution level
    :param idx: Index of the pixel to decompose
    :param mvec: SHD coefficients of the TF bin to be analysed
    :param k: Wave number (2*pi*f/c)
    :param ra: Radius of the microphone array
    :param Nmax: Maximum order of the SHD to be used
    :param pixeldict: Precalculated eigenvalue and eigenvector matrices (see loadpixbasis in utils.py)
    :return: New nodes after expanding the representation (if new entropy is less than the old entropy)
    """
    tn = tr.copy()
    newnodes = defaultdict()
    try:
        tn.pop((level, idx))
        lv = level + 1
        ids = 4 * idx
        A = hp.nside2pixarea(2 ** lv)
        for ind in range(4):
            w, V = pixeldict[(lv, ids + ind)]
            SRP = srpd(mvec, k, ra, Nmax, w, V)
            tn[(lv, ids + ind)] = SRP / A
        en = treeEntropy(tn)
        if en < er:
            cld = children(level, idx)
            for ind in range(4):
                newnodes[cld[ind]] = tn[cld[ind]]
    except KeyError:
        newnodes = defaultdict()
    return newnodes


def preprocessinput(insig, Nmax, NFFT, olap):
    """
    Preprocess the input to represent in TF domain and then obtain the SHD

    :param insig: (32 x length) matrix containing the audio signal channels from em32
    :param Nmax: Maximum SHD order
    :param NFFT: FFT size
    :param olap: (NFFT / olap) is the hop_size of the STFT
    :return: SHD of the STFT and the STFT of the microphone array recordings

    Notes:
    1. Currently accepts 32 channels (Eigenmike em32) is hard coded into this method.
    2. Uses STFT from madmom library since onset detection is also done by it. This can be replaced with any other TF
    bin selection method (e.g. direct-path dominance, RENT, soundfield directivity etc.).
    """
    P = []
    for ind in tqdm.trange(32, desc="Preprocessing channels"):
        ss = mm.audio.stft.STFT(insig[ind, :], frame_size=NFFT, hop_size=NFFT / olap, fft_size=NFFT)
        P.append(np.array(ss))
    Pnm = shd_all(P, Nmax)  # Calculate ths SHD from the STFTs
    return Pnm, P


def binstoprocess_onset(Pnm0, sg, fL, fH, Fs, NFFT, olap):
    """
    Select the bins with energy higher than RMS mean

    :param Pnm0: 0th-order SHD coefficient STFT matrix
    :param sg: 32xlength audio signals from em32
    :param fL: Lower frequency bound
    :param fH: Higher frequency bound
    :param Fs:  Sampling rate (Hz)
    :param NFFT: FFT size
    :param olap: (NFFT / olap) is the hop_size of the STFT
    :return: Time and frequency indices of the bins that are selected for HiGRID DOA estimation
    """

    em32 = EigenmikeEM32()
    estr = em32.returnAsStruct()
    wts = estr['weights']
    f1 = int(fL / Fs * NFFT)
    f2 = int(fH / Fs * NFFT)
    sgo = np.zeros(sg[0].shape) * 1j
    for ind in range(32):  # Calculate omnidirectional response
        sgo += sg[ind] * wts[ind]

    ls = mm.audio.spectrogram.LogarithmicFilteredSpectrogram(sgo, frame_size=NFFT, hop_size=NFFT / olap, fft_size=NFFT,
                                                             num_bands=24, sample_rate=Fs)
    os = mm.features.onsets.superflux(ls)  # These are the onsets
    os = os / max(os)  # These are the onsets in the analysed recording
    THR = np.sqrt(np.mean(np.array(os) ** 2))
    pk = indexes(os, thres=THR, min_dist=2)
    idx = []
    idy = []
    for ind in pk:
        vect = np.abs(np.abs(np.array(Pnm0[ind, f1:f2 + 1])))
        vb = baseline(vect, deg=2)
        vect, _ = scale(vect - vb)
        THRv = 0.1
        ids = indexes(vect, thres=THRv, min_dist=1)
        for knd in ids:
            idx.append(ind)
            idy.append(f1 + knd)
    return idx, idy


def binstoprocess_dpd(P, fL, fH, Fs, NFFT, Jnu, Jtau, Ndec, thr):
    """
    Select the bins with unit effective rank of the spatial autocorrelation matrix

    :param P: STFTs of the RSMA (e.g. em32) recordings
    :param fL: Lower frequency bound (Hz)
    :param fH: Upper frequency bound (Hz)
    :param Fs: Sampling rate (Hz)
    :param NFFT: FFT size
    :param Jnu: Frequency smoothing window size
    :param Jtau: Time smoothing window size
    :param Ndec: SHD order
    :param thr: DPD threshold
    :return: Time and frequency indices of the bins that are selected for HiGRID DOA estimation
    """

    em32 = EigenmikeEM32().returnAsStruct()
    fimin = int(round(fL / Fs * NFFT))
    fimax = int(round(fH / Fs * NFFT))

    # W, Y = getWY(em32, Ndec)

    Bmat = getBmat(em32, fimin, fimax + Jnu, NFFT, Fs, Ndec)
    Anm = getAnm(P, em32, Bmat, fimin, fimax + Jnu, Ndec)

    imax = Anm[0].shape[0]

    idx = []
    idy = []
    for tind in tqdm.trange(imax - Jtau, desc="DPD bin selection     "):
        for find in range(fimin, fimax):
            flag = dpd(Anm, Ndec, find, tind, Jtau, Jnu, thr)
            if flag:
                idx.append(tind)
                idy.append(find)
    return idx, idy


def higrid(mvec, freq, Ndec, tLevel, pixeldict):
    """
    Return the multiresolution HiGRID decomposition for a single TF bin

    :param mvec: SHD coefficients of the TF bin to be analysed
    :param freq: Frequency (in Hz)
    :param Ndec: SHD decomposition level
    :param tLevel: Resolution level
    :param pixeldict: Precalculated eigenvalue and eigenvector matrices (see loadpixbasis in utils.py)
    :return: Multiresolution HiGRID decomposition for a single TF bin
    """
    trr = treeroott(mvec, 2 * np.pi * freq / 340, 4.2e-2, Ndec, pixeldict)
    a = defaultdict()

    for lvl in range(0, tLevel):
        a.clear()
        Nidx = hp.nside2npix(2 ** lvl)
        er = treeEntropy(trr)
        bidx = range(Nidx)
        shuffle(list(bidx))
        for idx in bidx:
            newnodes = selecttreet(trr, er, lvl, idx, mvec, 2 * np.pi * freq / 340, 4.2e-2, Ndec, pixeldict)
            a.update(newnodes)
        for node in a:
            levc = node[0]
            idxc = node[1]
            par = parent(levc, idxc)
            if par in trr:
                trr.pop(par)
            trr[node] = a[node]
    return trr


def higridestimate(sg, pixeldict, maxnum=1000, Fs=48000., Ndec=4, NFFT=1024, olap=16, tLevel=3, dpdflag=False, fL=2608.,
                   fH=5216., thr=6):
    """
    Return the DOA estimates for multiple sources using HiGRID

    :param dpdflag:
    :param sg: Individual audio channels (32 channels) of em32
    :param pixeldict:  Precalculated eigenvalue and eigenvector matrices (see loadpixbasis in utils.py)
    :param maxnum: Maximum number of bins to be processed (default = 1000)
    :param Fs: Sampling frequency in Hz (default = 48000)
    :param Ndec: SHD decomposition order (default = 4)
    :param NFFT: FFT size (default = 1024)
    :param olap: (NFFT / olap) is the hop_size of the STFT (default = 16)
    :param tLevel: Tree decomposition level (default = 3)
    :param dpdflag: If True use direct-path dominance to select TF bins (default = False)
    :param fL: Lower frequency bound (in Hz) (default = 2608.)
    :param fH: Upper frequency bound (in Hz) (default = 5216.)
    :param thr: DPD threshold (default = 6)
    :return: Azimuth and elevation angles of the identified sources with respect to the array front direction (in radians)
    """

    Pnm, P = preprocessinput(sg, Ndec, NFFT, olap)

    if dpdflag:
        idx, idy = binstoprocess_dpd(P, fL, fH, Fs, NFFT, 25, 4, 3, thr)
    else:
        idx, idy = binstoprocess_onset(Pnm[0], sg, fL, fH, Fs, NFFT, olap)

    idx, idy = selectsome(idx, idy, maxnum)
    numbinsp = len(idx)

    vects = dict(); vects['th'] = []; vects['ph'] = []

    for ind in tqdm.trange(numbinsp, desc="HiGRID DOA Estimation "):

        mvec, NDec = selectbinindx(idy[ind], idx[ind], Pnm)
        freq = float(idy[ind]) * Fs / NFFT
        tr = higrid(mvec, freq, NDec, tLevel, pixeldict)
        cls = tc.treeCluster(tr, tLevel)
        cls.pop('type')

        for knd in cls.keys():
            nv = np.zeros(3)
            for item in cls[knd]:
                nv += node2vec(item[0], item[1])
            r, th, ph = cart2sph(nv[0], nv[1], nv[2])
            if ph < 0:
                ph = ph + 2 * np.pi
            vects['th'].append(th)
            vects['ph'].append(ph)

    thedge = np.linspace(0, np.pi, 181)
    phedge = np.linspace(0, 2 * np.pi, 361)[0:360]
    H, xe, ye = np.histogram2d(vects['th'], vects['ph'], bins=(thedge, phedge))

    Hn = H.T

    Hn[np.where(Hn < 2)] = 0
    Hd = median_filter(Hn, 3)

    Hg = gaussian_filter(Hd, sigma=3)

    trfinal = histtotree(Hg, thedge, phedge, tLevel)
    cls = tc.treeCluster(trfinal, tLevel)
    cls.pop('type')
    phloc = []
    thloc = []

    for knd in cls.keys():
        nv = np.zeros(3)
        for item in cls[knd]:
            nv += node2vec(item[0], item[1])
        r, th, ph = cart2sph(nv[0], nv[1], nv[2])
        if ph < 0:
            ph = ph + 2 * np.pi
        phloc.append(ph)
        thloc.append(th)

    return np.array(thloc), np.array(phloc)


if __name__ == '__main__':

    print("Please see higrid_example.py for an example usage")



