import struct
import pickle as pkl
from collections import defaultdict
from os import getcwd
import healpy as hp
import numpy as np
import wave
from scipy import signal as sp, special as sp


def wavread(wave_file):
    """
    Returns the contents of a wave file

    :param wave_file: Path to the wave_file to be read
    :return: (signal, sampling rate, number of channels)

    NOTE: Wavread solution was adapted from https://bit.ly/2Ubs9Jp
    """

    w = wave.open(wave_file)
    astr = w.readframes(w.getnframes())
    nchan = w.getnchannels()
    totsm = w.getnframes()
    sig = np.zeros((nchan, totsm))
    a = struct.unpack("%ih" % (w.getnframes() * w.getnchannels()), astr)
    a = [float(val) / pow(2, 15) for val in a]
    for ind in range(nchan):
        b = a[ind::nchan]
        sig[ind] = b

    sig = np.transpose(sig)
    fs = w.getframerate()
    w.close()
    return sig, fs, nchan


def processirs(irs, monosnd):
    """
    Returns an emulated recording of em32 using acoustic impulse responses and an anechoic sound signal

    :param irs: Acoustic impulse responses obtained using em32
    :param monosnd: Monophonic sound signal to be convolved with the AIRs
    :return: 32-channels emulated recording of em32 as a numpy array
    """
    nmirs = len(irs)
    lnsnd = len(monosnd)
    lnirs = np.size(irs) / nmirs

    out = np.zeros((nmirs, lnsnd + lnirs - 1))

    for ind in range(nmirs):
        ir = irs[ind].flatten()
        snd = monosnd.flatten()
        ot = sp.fftconvolve(ir, snd)
        out[ind] = ot
    return out


def node2vec(level, idx):
    """
    Converts the centre direction of a Healpix pixel to a unit vector

    :param level: Resolution level
    :param idx: Index of the pixel
    :return: 3x1 array containing the unit vector
    """
    th, ph = hp.pix2ang(2 ** level, idx, nest=True)
    # print(th/np.pi*180, ph/np.pi*180)
    x = np.sin(th) * np.cos(ph)
    y = np.sin(th) * np.sin(ph)
    z = np.cos(th)
    return np.array([x, y, z])


def selectsome(idx, idy, maxnum):
    """
    Randomly selects a given number of time-frequency bins

    :param idx: List containing frequency indices of selected bins
    :param idy: List containing time indices of selected bins
    :param maxnum: Maximum number of bins to select (in None, return all indices)
    :return:
    """
    if maxnum == None or len(idx) < maxnum:
        return idx, idy
    else:
        ids = np.random.randint(0, len(idx)-1, maxnum)
        idx = list(np.array(idx)[ids])
        idy = list(np.array(idy)[ids])
        return idx, idy


def cart2sph(x, y, z):
    """
    r, th, ph = cart2sph(x, y, z)

    Return the spherical coordinate representation of point(s) given in
    Cartesian coordinates

    As usual r is the radius, th is the elevation angle defined from the
    positive z axis and ph is the azimuth angle defined from the positive x axis
    """
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    th = np.arccos(z / r)
    ph = np.arctan2(y, x)
    return r, th, ph


def sph2cart(r, th, ph):
    """
    Converts vector in spherical coordinates to Cartesian coordinates

    :param r: Radius
    :param th: Azimuth angle
    :param ph: Inclination angle
    :return: Vector in Cartesian coordinates
    """
    x = r * np.sin(th) * np.cos(ph)
    y = r * np.sin(th) * np.sin(ph)
    z = r * np.cos(th)
    u = np.array([x, y, z])
    return u


def selectbinindx(fidx, tidx, Pnm, Ndec = 4):
    """
    Returns the SHD vectors of a given time-frequency bin

    :param fidx: Frequency index
    :param tidx: Time index
    :param Pnm: List of SHD-STFT matrices
    :param Ndec: SHD order (default = 4)
    :return: SHD vector, SHD order
    """
    mvec = np.zeros((Ndec + 1) ** 2) * 1j
    for ind in range((Ndec + 1) ** 2):
        M = Pnm[ind]
        mvec[ind] = M[tidx, fidx]
    return mvec, Ndec

def loadpixbasis():
    """
    Returns the

    :return:
    """
    fl = open(getcwd() + '/data/cdata/pixelbasis_red.pkl', 'rb')
    pd = pkl.load(fl, encoding='latin1')
    fl.close()
    return pd


def fulltree(tLevel=3, val=0):
    """
    Create a defaultdict containing all the pixels in a healpix grid, initialised with the default value

    :param tLevel: Healpix resolution level (default = 3)
    :param val: Value of each pixel (default = 0)
    :return: defaultdict containing all the pixels in the healpix grid at the given resolution
    """
    tree = defaultdict()
    nside = hp.nside2npix(2**tLevel)
    for ind in range(nside):
        tree[(tLevel, ind)] = val
    return tree


def histtotree(H, the, phe, tLevel):
    """
    Convert 2D DOA histogram to a healpix tree representation

    :param H: 2D DOA histogram matrix
    :param the: Array of azimuth angles corresponding to columns of the 2D histogram matrix
    :param phe: Array of inclination angles corresponding to columns of the 2D histogram matrix
    :param tLevel: Healpix resolution level (default = 3)
    :return: defaultdict containing healpix grid st the given resolution containing the 2D histogram
    """
    ft = defaultdict()
    for ind in range(len(the)-1):
        for jnd in range(len(phe)-1):
            th = (the[ind]+the[ind+1]) / 2
            ph = (phe[jnd]+phe[jnd+1]) / 2
            idpix = hp.ang2pix(2**tLevel, th, ph, nest=True)
            if (tLevel, idpix) in ft:
                ft[(tLevel, idpix)] += H[jnd, ind]
            else:
                ft[(tLevel, idpix)] = H[jnd, ind]
    return ft

def sph_jnyn(N, kr):
    '''
    Returns spherical Bessel functions of the first (jn) and second kind (yn) and their derivatives

    :param N: Function order
    :param kr: Argument
    :return: jn, jn', yn, yn'

    NOTE: Emulates the behaviour of sph_jnyn() in early versions of scipy (< 1.0.0).

    '''
    jn = np.zeros(N+1)
    jnp = np.zeros(N+1)
    yn = np.zeros(N+1)
    ynp = np.zeros(N+1)
    for n in range(N+1):
        jn[n] = sp.spherical_jn(n, kr)
        jnp[n] = sp.spherical_jn(n, kr, derivative=True)
        yn[n] = sp.spherical_yn(n, kr)
        ynp[n] = sp.spherical_yn(n, kr, derivative=True)
    return jn, jnp, yn, ynp
