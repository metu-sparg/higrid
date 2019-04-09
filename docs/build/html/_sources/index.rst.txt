.. higrid documentation master file, created by
   sphinx-quickstart on Tue Apr  9 12:11:02 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

====================
HiGRID documentation
====================

`higrid` is a Python package that implements the hierarchical grid refinement (HiGRID) direction-of-arrival estimation algorithm for rigid spherical microphone arrays. The algorithm is based on the calculation of steered response power density (SRPD) maps, and spatial entropy-based (multiple) peak detection.

The method was developed by the researchers in METU Spatial Audio Research Group http://www.sparglab.org. The technical details are available in the following papers. Please cite our papers if you want to use the code and/or the data provided in this package in your research.

*Coteli, M. B., Olgun, O., and Hacihabiboglu, H. (2018). Multiple Sound Source Localization With Steered Response Power Density and Hierarchical Grid Refinement. IEEE/ACM Trans. Audio, Speech and Language Process., 26(11), 2215-2229.* https://ieeexplore.ieee.org/document/8418732

*Olgun, O. and Hacihabiboglu, H., (2018) "Localization of Multiple Sources in the Spherical Harmonic Domain with Hierarchical Grid Refinement and EB-MUSIC". In 2018 16th Int. Workshop on Acoust. Signal Enhancement (IWAENC-18) (pp. 101-105), Tokyo, Japan.* https://ieeexplore.ieee.org/abstract/document/8521365 ; also available at https://bit.ly/2I81eru

Installation
============

Due to one of the key external packages (i.e. healpy) working on Linux and MacOS only (and not Windows), higrid will also only work on Linuz and MacOS, only. The current version was tested on MacOS 10.14.1 with Python 3.6.1.

Prerequisites
-------------

In order to install and use the `higrid` you need Python 3.x on your system. A few tweaks would be necessary to make it work with Python 2.7 (incompatibilities are due to scipy version, as well as default parameters for pickle implementation in 2.7 and defaultdict member functions in 2.7). The following packages are also needed:

* numpy http://www.numpy.org
* scipy https://www.scipy.org
* healpy https://github.com/healpy/healpy
* PeakUtils https://bitbucket.org/lucashnegri/peakutils
* tqdm https://www.github.com/tqdm/tqdm
* madmom https://github.com/CPJKU/madmom

Please also check the requirements.txt file at https://github.com/metu-sparg/higrid/blob/master/requirements.txt for the specific version requirements.

Installation from `pip3` Package
--------------------------------

Using `pip3` is the easiest way to install `higrid` :

.. code-block:: python

   pip3 -install higrid

Example usage
=============

**Please see the file `higrid_example.py` for a more complete example.**

.. code-block:: python
   :linenos:

   import higrid as hg

   testinstance = set([(1, 1, 2), (5, 1, 2), (5, 5, 2), (1, 5, 2)])
   sg = hg.composescene(
        ['music/mahler_vl1a_6.wav', 'music/mahler_vl1b_6.wav', 'music/mahler_vl2a_6.wav', 'music/mahler_vl2b_6.wav'],
        testinstance, (0, 192000))
   pd = hg.loadpixbasis()
   th, ph = hg.higridestimate(sg, pd, maxnum=1000, Fs=48000, Ndec=4, NFFT=1024, olap=16, tLevel=3, dpdflag=True,
                               fL=2608., fH=5216., thr=6)

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   modules

==================
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
