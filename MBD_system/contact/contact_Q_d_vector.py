'''
Created on 25. apr. 2015

@author: lskrinjar
'''
import itertools
import numpy as np

class Contact_Q_d_vector(object):
    '''
    classdocs
    '''
    __id = itertools.cycle([0, 1])

    def __init__(self, C_qq, C_ss, C_qs, dq, ds):
        '''
        Constructor
        '''
        self.id = self.__id.next()

        self.vector = 2 * np.dot(C_qs, dq) * ds + np.dot(C_qq, dq ** 2) + np.dot(C_ss, ds ** 2)  # 
