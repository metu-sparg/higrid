import math

def parent(level, idx):
    """
    Return the parent of a given healpix pixel (in nested format)

    :param level: Resolution level
    :param idx: Index at the given level
    :return: Tuple (lvl, idp) including the level and the index of the parent pixel
    """
    assert idx < 12 * 2 ** (2 * level)
    lvl = level - 1
    idp = int(math.floor(idx / 4))
    return (lvl, idp)


def parents(level, idx):
    """
    Return all the (grand-)parents of the Healpix pixel idx at level (in nested format)

    :param level: Resolution level
    :param idx: Pixel index
    :return: All the parents of the pixel
    """
    assert idx < 12 * 2 ** (2 * level)
    plpairs = []
    for ind in range(level, 0, -1):
        idx = int(math.floor(idx / 4))
        plpairs.append(tuple((ind - 1, idx)))
        level -= 1
    return plpairs[::-1]


def children(level, idx):
    """
    Return all the children of the Healpix pixel idx at level (in nested format)

    :param level: Resolution level
    :param idx: Pixel index
    :return: All the parents of the pixel
    """
    chld = []
    for ind in range(4):
        chld.append((level + 1, 4 * idx + ind))
    return chld


def siblings(level, idx):
    """
    Return the siblings of the Healpix pixel idx at level (in nested format)

    :param level: Resolution level
    :param idx: Pixel index
    :return: Return the siblings of the pixel
    """
    sibl = children(level - 1, idx / 4)
    return sibl