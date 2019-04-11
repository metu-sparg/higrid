# higrid

[![DOI](https://zenodo.org/badge/180388725.svg)](https://zenodo.org/badge/latestdoi/180388725)

`higrid` is a Python package that implements the hierarchical grid refinement (HiGRID) direction-of-arrival estimation algorithm based on steered response power density (SRPD) maps, and spatial entropy-based (multiple) peak detection.

The method was developed by the researchers in METU Spatial Audio Research Group ([http://www.sparglab.org](http://www.sparglab.org)). The technical details are available in the following papers. Please cite our papers if you want to use the code and/or the data provided in this package.

*Coteli, M. B., Olgun, O., and Hacihabiboglu, H. (2018). Multiple Sound Source Localization With Steered Response Power Density and Hierarchical Grid Refinement. IEEE/ACM Trans. Audio, Speech and Language Process., 26(11), 2215-2229.* [[Link]](https://ieeexplore.ieee.org/document/8418732)

*Olgun, O. and Hacihabiboglu, H., (2018) "Localization of Multiple Sources in the Spherical Harmonic Domain with Hierarchical Grid Refinement and EB-MUSIC". In 2018 16th Int. Workshop on Acoust. Signal Enhancement (IWAENC-18) (pp. 101-105), Tokyo, Japan.* [[Link 1]](https://ieeexplore.ieee.org/abstract/document/8521365), [[Link 2]](https://www.researchgate.net/publication/328137952_Localization_of_Multiple_Sources_In_The_Spherical_Harmonic_Domain_With_Hierarchical_Grid_Refinement_and_EB-MUSIC?_sg=qSld-Z7cru_p-KEcz20GJOLuGgy97dycxdQ2aCPWZUP87181q9OmpUU1U_0uOmxo8TAkDZPcpN23BwajIj1GxM6Irm4)


## Documentation

Documentation of the package is available online [http://higrid.readthedocs.org](http://higrid.readthedocs.org)

## License

The package has different licenses for the source code and for the data files.

### Source Code

Unless indicated otherwise, all source code files are published under the BSD license. For details, please see the LICENSE.txt file.

### Acoustic Impulse Responses (AIR)

Unless indicated otherwise, acoustic impulse response (AIR) data files are distributed under the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 license](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode). This dataset along with a documentation is also provided on Zenodo (https://doi.org/10.5281/zenodo.2635758)

### Audio recordings

Three sets of audio recordings are provided with different licenses for each.

#### Anechoic Speech Signals

Anechoic speech recordings were made at METU SPARG Audio Lab by [Ms Özgen Demirkaplan](https://www.researchgate.net/profile/Oezgen_Demirkaplan) and used in the following paper:

*Demirkaplan, O., and Hacihabiboglu, H., (2019) The effect of interpersonal familiarity on auditory distance perception of reverberant speech, submitted to J. Acoust. Soc. Am. (revised version under review)*

These recordings are distributed under the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 license](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode).

#### Anechoic Music Signals

These anechoic orchestra recordings are used to demonstrate the utiliy of the HiGRID model for near-coherent sources. These are excerpts from anechoic recordings kindly provided by Jukka Pätynen, Ville Pulkki, and Tapio Lokki from Aalto University. 

*Pätynen, J., Pulkki, V., and Lokki, T., "Anechoic recording system for symphony orchestra," Acta Acustica united with Acustica, vol. 94, nr. 6, pp. 856-865, November/December 2008.*[http://dx.doi.org/10.3813/AAA.918104](http://dx.doi.org/10.3813/AAA.918104)
 
For the intellectual rights and the distribution policy of the audio recordings in this dataset contact Aalto University. For more information about the original anechoic recordings we refer to the [web page](https://users.aalto.fi/~ktlokki/Sinfrec/sinfrec.html) and the associated publication given above.

#### Real Music Recordings with Eigenmike em32

These recordings were made during the rehearsals of [Nemeth Quartet](http://www.nemethquartet.com) at the [Erimtan Museum](http://erimtanmuseum.org/) recital hall in Ankara, Turkey on March, 7. 2017 and used here to demonstrate the capabilities of HiGRID with a real recording. **The recordings are copyrighted and should not used without explicit permission.**

## Installation

### Prerequisites
In order to install and use the `higrid` you need Python 3.x on your system. A few tweaks would be necessary to make it work with Python 2.7 (incompatibilities are due to scipy version, as well as default parameters for pickle implementation in 2.7 and defaultdict member functions in 2.7). The following packages are also needed:

* [`numpy`](http://www.numpy.org)
* [`scipy`](https://www.scipy.org)
* [`healpy`](https://github.com/healpy/healpy)
* [`PeakUtils`](https://bitbucket.org/lucashnegri/peakutils)
* [`tqdm`](https://www.github.com/tqdm/tqdm)
* [`madmom`](https://github.com/CPJKU/madmom)

Please also check the [requirements.txt](https://github.com/metu-sparg/higrid/blob/master/requirements.txt) file for the specific version requirements.

### Installation using `pip3`

Using `pip3` is the easiest way to install: `pip3 install higrid`

### Installation from source

If you want to download and use the source code from Github:

`git clone --recursive https://github.com/metu-sparg/higrid.git`

## Acknowledgments

The code and data provided here was created in a research project supported by the [Turkish Scientific and Technological Research Council](http://www.tubitak.gov.tr) (TUBITAK) via Research Grant 113E513 “Spatial Audio Reproduction Using Analysis-based Synthesis Methods” (2014-2018).
