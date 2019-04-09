from os import getcwd, listdir
import numpy as np
from scipy import signal as sp

from higrid.utils import wavread

def emulatescene(insig, gain, irspath):
    """
    Emulates a scene by convolving a source input signal with em32 AIRs

    :param insig: (Single-channel) input signal
    :param gain: Gain (scalar)
    :param irspath: Path to the AIRs to be used
    :return: 32 channels of audio from an emulated em32 recording.
    """
    dr = listdir(irspath)
    dr = sorted(dr)
    wv = wavread(irspath + '/' + dr[0])
    ir = np.zeros((32,wv[0].shape[0]))
    for ind in range(32):
        wv = wavread(irspath + '/' + dr[ind])
        ir[ind,:] += wv[0].reshape((wv[0].shape[0]))

    sz = len(insig)
    out = np.zeros((32, sz))
    for ind in range(32):
        out[ind,:] = sp.fftconvolve(gain * insig, ir[ind,:], mode='same')
    return out

def emptyscene(sha = (32, 48000)):
    """
    Returns an empty scene

    :param sha: Tuple containing number of channels and number of samples (nchan, samples) (default = (32, 48000))
    :return: An empty scene containing nchan channels and the given number of samples (empty numpy array)
    """
    sge = np.zeros(sha)
    return sge

def combinescene(sg1, sg2):
    """
    Linearly combines two scenes pertaining to em32 recordings

    :param sg1: Scene 1 (32 x samples numpy array)
    :param sg2: Scene 2 (32 x samples numpy array)
    :return: Combined scene (32 x samples numpy array)
    """
    sgo = np.zeros(sg1.shape)
    for ind in range(32):
        sgo[ind,:] = sg1[ind,:] + sg2[ind,:]
    return sgo

def composescene(filelist, dirset, samples=(0, 96000), roomstr='ii-s05'):
    """
    Compose an emulated scene using a number of anechoic sound signals and measured AIRs

    :param filelist: List of files to be used
    :param dirset: Set containing tuples with (X, Y, Z) as the AIR indices
    :param samples: Start and end points of samples to be prococessed as a tuple (sstart, send)
    :param roomstr: Used to select from a specific directory (default is 'ii-s05' as we only provided AIRs for that room)
    :return: 32 channels of audio from an emulated em32 recording.
    """
    drset = dirset.copy()
    assert len(filelist) == len(drset)
    numsamp = samples[1] - samples[0]
    sgo = emptyscene((32, numsamp))
    for item in filelist:
        dr = drset.pop()
        drtxt = str(dr[0]) + str(dr[1]) + str(dr[2])
        snd = wavread(getcwd() + '/data/sdata/anechoic/' + item)
        snd = snd[0].reshape((snd[0].shape[0]))[samples[0]:samples[1]]
        gain = np.sqrt((dr[0]-3.0)**2 + (dr[1]-3.0)**2 + ((dr[2]-2.0)*0.6)**2)
        sg = emulatescene(snd, gain, getcwd() +'/data/rdata/'+ roomstr + '/' + drtxt)
        sgo = combinescene(sgo, sg)
    return sgo

def realrec(dirpath, prefix, samples):
    """
    Create a scene from real em32 recordings

    :param dirpath: Path containing the em32 recordings
    :param prefix:
    :param samples: Number of samples to use
    :return: 32 channel em32 recording
    """
    sg = emptyscene((32,samples))
    for ind in range(32):
        snd = wavread(dirpath + prefix + str(ind+1) + '.wav')
        snd = snd[0].reshape((snd[0].shape[0]))[0:samples]
        sg[ind,:] = snd
    return sg