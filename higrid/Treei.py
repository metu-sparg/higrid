from collections import defaultdict
import healpy as hp

from higridestimate import srpd
from treeutils import children


def selecttreet(tr, er, level, idx, mvec, k, ra, Nmax, pixeldict):
    """
    Expand a given node based on its effect on spatial entropy

    :param tr: Healpix representation of the SRPD map
    :param er: Initial spatial entropy of the tree
    :param level: Resolution level
    :param idx: Pixel index
    :param mvec: SHD of of the STFT bin
    :param k: Wave number
    :param ra: Array radius
    :param Nmax: SHD order
    :param pixeldict: Precalculated eigenvalue and eigenvector matrices (see loadpixbasis in utils.py)
    :return:
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