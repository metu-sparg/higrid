import numpy as np
from os import getcwd

import higrid as hg



"""
higrid_example.py: Presents an example usage for HiGRID DOA estimation. Please see the comments in the code. 

NOTES:  1.  The code was tested on Mac OS 10.14.1 with Python 3.6. 
        2.  A known problem with other OSs (e.g. Ubuntu) mainly happens because of the differences in the directory 
            structure, which effectively results in problems with channel ordering.
        2.  HiGRID relies heavily on the healpy package which is only supported on Linux and Mac OS.
        3.  Please make sure that the individual microphone array channels are read in proper order.
"""

if __name__=='__main__':

    """
       1a.   Select a number of source directions from the measurement grid
             The example below is for positions ar 45, 135, 225, 315 degrees azimuth in the horizontal plane
    """

    testinstance = set([(1, 1, 2), (5, 1, 2), (5, 5, 2), (1, 5, 2)])

    """1b.   Uncomment lines below to create emulations for testing the algorithm using the measurement grid"""

    # mg = hg.measgrid()
    # testinstance = hg.selectset(mg, 4, 3, 8, np.pi / 4, 1)[0][0]

    """""""# 1c.   Compose an emulated scene..."""

    """i) ...using near-coherent music signals"""
    sg = hg.composescene(
        ['music/mahler_vl1a_6.wav', 'music/mahler_vl1b_6.wav', 'music/mahler_vl2a_6.wav', 'music/mahler_vl2b_6.wav'],
        testinstance, (0, 192000))

    """ii) OR comment above and uncomment below to use speech signals"""
    # sg = hg.composescene(
    #    ['speech/male1.wav', 'speech/female2.wav', 'speech/male3.wav', 'speech/female4.wav'],
    #    testinstance[0][0], (0, 192000))

    """1d.   Alternatively, comment everything above and uncomment the following to use real recordings."""
    # sg = hg.realrec(getcwd() + '/data/sdata/real/', 'Quartet', 192000)

    """2.    Load the pickle file that contains the precomputed eigenvalues and eigenvector matrices"""
    pd = hg.loadpixbasis()

    """
       3.    Calculate the DOA estimates using HiGRID
             NOTES:
             1.  Selection of DPD threshold is not yet implemented in this version of the package.
             2.  For real recordings see:
                 Coteli, Mert Burkay, Orhun Olgun, and Huseyin Hacihabiboglu. "Multiple Sound Source Localization With
                 Steered Response Power Density and Hierarchical Grid Refinement." IEEE/ACM Transactions on Audio, Speech
                 and Language Processing (TASLP) 26, no. 11 (2018): 2215-2229.
             3.  Data is used in this example have different license/copyright/usage rights. Please see the Github 
                 repository for details.
    """

    th, ph = hg.higridestimate(sg, pd, maxnum=1000, Fs=48000, Ndec=4, NFFT=1024, olap=16, tLevel=3, dpdflag=True,
                               fL=2608., fH=5216., thr=6)

    """4. Convert to degrees and show the results"""
    azims = ph / np.pi * 180.0
    elevs = th / np.pi * 180.0
    numsrc = len(azims)

    print(" ")
    print("HiGRID DOA ESTIMATES")
    print("--------------------")

    for jnd in range(numsrc):
        print("Source " + str(jnd + 1) + ": " + "{:3.2f}".format(azims[jnd]) + "   " + "{:3.2f}".format(elevs[jnd]))
