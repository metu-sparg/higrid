import healpy as hp
import numpy as np
from collections import defaultdict


def preprocess(trr, tLevel):
    """
    Eliminate leaf nodes that are below the mean

    :param trr: Healpix tree (as a defaultdict)
    :param tLevel: Level at which the leaf nodes are
    :return: Tree as a defaultdict after mean thresholding
    """
    tr = trr.copy()
    hlev = [item for item in tr if item[0] == tLevel] # Get the highest level leaf nodes
    travr = np.mean([abs(tr[item]) for item in tr.keys()])
    for item in hlev:
        if abs(tr[item]) < travr:
            tr.pop(item)
    return tr

def hpneighbouridx(level, idx):
    """
    Neighbours of a Healpix pixel (in nested format)

    :param level: Resolution level
    :param idx: Index of the Healpix pixel
    :return: Neighbouring nodes of the pixel
    """
    nside = 2**level
    th, ph = hp.pix2ang(nside, idx, nest=True)
    ninds = hp.get_all_neighbours(nside, th, ph, nest=True)
    nbrs = []
    lev = int(np.log2(nside))
    for item in ninds:
        nbrs.append((lev, item))
    return nbrs

def neighboursofset(cset):
    """
    Return all the neighbours of a set of Healpix pixels (in nested format)
    :param cset:
    :return:
    """
    cbar = set()
    for item in cset:
        nbrs = hpneighbouridx(item[0], item[1])
        cbar = cbar.union(set(nbrs))
    return cbar


def treeCluster(trr, tLevel):
    """
    Neighbouring Nodes Labelling (NNL) for clustering pixels

    :param trr: Healpix tree (as a defaultdict)
    :param tLevel: Level at which the leaf nodes are
    :return: Healpix pixel clusters that belong to a single source
    """
    tr = trr.copy()
    tr = preprocess(tr, tLevel)
    hlev = [item for item in tr if item[0] == tLevel] # Get the highest level leaf nodes
    clusters = defaultdict(type=set)
    cind = 1
    allset = set(hlev)

    while allset: # Continue the process until all the nodes are consumed
        flag = 1
        clusters[cind] = set()
        n = allset.pop()
        clusters[cind].add(n)
        while flag:
            cbar = neighboursofset(clusters[cind]) # Find all neighbours of the cluster
            ints = cbar.intersection(allset) # Find the neighbours which are in the tree
            if ints: # If there is at least one neighbour
                clusters[cind] = clusters[cind].union(ints)
                for item in ints:
                    allset.remove(item)
            else:
                cind += 1
                flag = 0
    return clusters