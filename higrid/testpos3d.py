import numpy as np
from collections import defaultdict
import random
import csv as csv
from os import getcwd


def vdist(a, b):
    """
    Returns the angle between two vectors

    :param a: Vector 1
    :param b: Vector 2
    :return: Normalised angle [-1,+1] between the two vectors
    """
    dst = np.arccos(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))) / np.pi
    return dst

def measgrid():
    """
    Returns the measurement positions and DRR values for emulations

    :return: Defaultdict containing tuple (xindex, yindex, zindex) as keys and DRR as values
    """
    fl = open(getcwd() + '/data/cdata/drr_3d.csv', 'rU')
    m = csv.reader(fl, delimiter=';')
    DRR = np.zeros((244,1))
    IDX = np.zeros((244,3))
    ind = 0
    for row in m:
        x = int(row[4])
        y = int(row[5])
        z = int(row[6])
        D = float(row[3])
        IDX[ind,:] = [y, x, z] #??
        DRR[ind] = D
        ind += 1

    sdict = defaultdict()
    mid = (3, 3, 2)
    for ind in range(244):
        sdict[(IDX[ind,0], IDX[ind,1], IDX[ind,2])] = DRR[ind]
    # Right above and right below the array. We don't use these positions.
    sdict.pop((3, 3, 0))
    sdict.pop((3, 3, 1))
    sdict.pop((3, 3, 3))
    sdict.pop((3, 3, 4))
    mgrid = defaultdict()
    mgrid['pos'] = sdict
    mgrid['origin'] = mid
    return mgrid


def angle(tupl1, tupl2, mtupl):
    """
    Returns the angle between two vectors

    :param tupl1: Tuple containing the coordinates of a point
    :param tupl2: Tuple containing the coordinates of another point
    :param mtupl: Tuple containing the coordinates of origin
    :return: Angle between the two vectors (in radians)
    """
    v1 = np.array(tupl1, dtype=float) - np.array(mtupl, dtype=float)
    v2 = np.array(tupl2, dtype=float) - np.array(mtupl, dtype=float)
    v1n = np.linalg.norm(v1)
    v2n = np.linalg.norm(v2)
    th = np.arccos(np.vdot(v1, v2) / (v1n * v2n + 1e-7))
    return th

def histclust(drrmat, dind, N, numthr):
    """
    Returns a number of clusters of sources positions based on their D/R ratios

    :param drrmat: Matrix containing coordinates of source oposition on the grid and its D/R ratio
    :param dind: Histogram bin to choose
    :param N: Number of bins to use
    :param numthr: At least this many sources in a given bin has to be present
    :return:
    """

    h, be = np.histogram(drrmat, N)
    assert dind < N
    assert h[dind] > numthr
    # print((be[1:]+be[0:-1])/2)
    # print('Bin size is: ' + str(be[1]-be[0]) + ' dB\n')
    ids = np.cumsum(h)
    id = np.concatenate(([0], ids))
    start = id[dind]
    stop = id[dind + 1]
    return start, stop

def select(grd, n=4, dind=3, dclust=8, thcond=np.pi/4):
    """
    Returns a single test scene instance

    :param grd: Measurement grid from measgrid()
    :param n: Number of sources to choose (default = 4)
    :param dind: Cluster index (default = 3)
    :param dclust: Number of bins to use in histogram for D/R ratio clustering (default = 8)
    :param thcond: Minimum angle between consecutive sources (in radians) (default = pi/4)
    :return: One test instance and the average D/R ratio at the selected positions
    """
    assert(dind < dclust)
    ddict = grd['pos'].copy()
    mid = grd['origin']
    smat = np.zeros((239, 4))
    ind = 0
    for item in ddict.keys():
        smat[ind,:] = [item[0], item[1], item[2], ddict[item][0]]
        ind += 1

    sind = np.argsort(smat[:,3])
    smat = smat[sind]
    start, stop = histclust(smat[:,3], dind, dclust, 4)
    smatn = smat[start:stop, :]
    drrmean = np.mean(smatn[:,3])

    sdict = defaultdict()
    for item in smatn:
        sdict[(item[0], item[1], item[2])] = ddict.pop((item[0], item[1], item[2]))

    src = []
    s0 = random.choice(sdict.keys()) # Select the first source
    s0 = (int(s0[0]), int(s0[1]), int(s0[2]))
    sdict.pop(s0)
    src.append(s0)

    while len(src)<n:
        sp = random.choice(sdict.keys())
        sp = (int(sp[0]), int(sp[1]), int(sp[2]))
        flag = True
        for item in src:
            ang = angle(item, sp, mid)
            flag = flag and (ang > thcond)
        if flag:
            src.append(sp)
    return src, drrmean



def selectset(grd, n=4, dind=3, dclust=8, thcond=np.pi/4, trialcount=1):
    """
    Returns a number of test scene instances with similar D/R ratio and a prescribed minimum angle between them

    :param grd: Measurement grid from measgrid()
    :param n: Number of sources to choose (default = 4)
    :param dind: Cluster index (default = 3)
    :param dclust: Number of bins to use in histogram for D/R ratio clustering (default = 8)
    :param thcond: Minimum angle between consecutive sources (in radians) (default = pi/4)
    :param trialcount: Number of randomly selected scenes (default = 1)
    :return: List containing trialcount number of test instances and the average D/R ratio at the selected positions
    """
    oset = []
    while len(oset) < trialcount:
        a, drrmean = select(grd, n, dind, dclust, thcond)
        a = set(a)
        flag = True
        for item in oset:
            s = set(a).difference(item)
            if len(s) == 0:
                flag = False
        if flag:
            oset.append(a)
    return oset, drrmean
