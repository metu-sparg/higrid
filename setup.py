from setuptools import setup

setup(
    name='higrid',
    version='0.1.0',
    packages=['higrid'],
    url='https://github.com/metu-sparg/higrid',
    license='BSD 3-Clause License',
    author='Huseyin Hacihabiboglu (METU SPARG)',
    author_email='husshho+pypi@gmail.com',
    description='Hierarchical Grid Refinement (HiGRID): Multiple source DOA Estimation with Spherical Microphone Arrays',
    keywords=['Localization', 'DOA estimation', 'Eigenmike', 'Acoustics', 'Microphone array', 'Sound'],
    install_requires=[  # I get to this in a second
        'numpy',
        'scipy',
        'madmom',
        'PeakUtils',
        'tqdm',
        'healpy',
        ] ,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3 :: Only'
        'Environment :: Console',
        'Operating System :: MacOS',
        'Operating System :: POSIX :: Linux',
        'Topic :: Multimedia :: Sound/Audio',
        'Topic :: Multimedia :: Sound/Audio :: Analysis',
        'Topic :: Multimedia :: Sound/Audio :: Capture/Recording',
    ]
    download_url='https://github.com/metu-sparg/higrid/archive/higrid_v01.tar.gz'
)
