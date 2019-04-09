import numpy as np

'''
Microphone.py: Contains the Microphone, MicrophoneArray and EigenmikeEM32 classes.
'''

class Microphone(object):
    '''
    Microphone class
    '''

    def __init__(self, name='Generic', version='1.0', direct='Omnidirectional'):
        '''
        Constructor

        :param name: Name of the microphone
        :param version: Version of the microphone
        :param direct: Directivity of the microphone (str)
        '''
        self._micname = name
        self._ver = version
        self._directivity = direct

    def getname(self):
        '''
        Getter for the name

        :return: Name (str) of the microphone object
        '''
        return self._micname

    def getversion(self):
        '''
        Getter for the version

        :return: Version (str) of the microphone object
        '''
        print(self._ver)

    def setname(self, name):
        '''
        Setter for the name attribute

        :param name: Name of the microphone (str)
        '''
        self._micname = name

    def setversion(self, version):
        '''
        Setter for the version attribute

        :param name: Version of the microphone (str)
        '''
        self._ver = version


class MicrophoneArray(Microphone):
    '''
    MicrophoneArray Class that inherits from the Microphone class
    '''

    def __init__(self, name, typ, version, direct):
        '''
        Constructor

        :param name: Name of the array
        :param typ: Type of the array (e.g. 'RSMA', 'OSMA', 'Linear')
        :param version: Version of the array (e.g. '1.0a')
        :param direct: Directivity of components (e.g. 'Onmidirectional')
        '''
        super(MicrophoneArray, self).__init__(name, version, direct)
        self._arraytype = None
        self.__arraytype = typ

    def gettype(self):
        '''
        Getter for array type

        :return: Type of the array (str)
        '''
        return self.__arraytype

    def settype(self, typ):
        '''
        Setter for array type

        :param typ: Type of the array (str)
        '''
        self.__arraytype = typ


class EigenmikeEM32(MicrophoneArray):
    '''
    Eigenmike em32 class that inherits from the MicrophoneArray class.
    '''
    def __init__(self):
        super(EigenmikeEM32, self).__init__('Eigenmike 32', 'Rigid Spherical', 17.0, 'Omni')
        self._numelements = 32

        self._thetas = np.array([69.0, 90.0, 111.0, 90.0, 32.0, 55.0,
                               90.0, 125.0, 148.0, 125.0, 90.0, 55.0, 21.0, 58.0,
                               121.0, 159.0, 69.0, 90.0, 111.0, 90.0, 32.0, 55.0,
                               90.0, 125.0, 148.0, 125.0, 90.0, 55.0, 21.0, 58.0,
                               122.0, 159.0]) / 180.0 * np.pi

        self._phis = np.array([0.0, 32.0, 0.0, 328.0, 0.0, 45.0, 69.0, 45.0, 0.0, 315.0,
                             291.0, 315.0, 91.0, 90.0, 90.0, 89.0, 180.0, 212.0, 180.0, 148.0, 180.0,
                             225.0, 249.0, 225.0, 180.0, 135.0, 111.0, 135.0, 269.0, 270.0, 270.0,
                             271.0]) / 180.0 * np.pi

        self._radius = 4.2e-2

        self._weights = np.ones(32)

        self._info = 'Eigenmike em32 needs to be calibrated using the software tool provided mh Acoustics before use.'

    def returnAsStruct(self):
        '''
        Returns the attributes of the Eigenmike em32 as a struct

        :return: dict object with the name, type, thetas, phis, radius, weights, version, numelements, directivity, info fields
        '''
        em32 = {'name': self._micname,
                'type': self._arraytype,
                'thetas': self._thetas,
                'phis': self._phis,
                'radius': self._radius,
                'weights': self._weights,
                'version': self._ver,
                'numelements': self._numelements,
                'directivity': self._directivity,
                'info': self._info}
        return em32



