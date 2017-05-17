"""
Created on 13. mar. 2014

@author: luka.skrinjar
"""
import logging
from pprint import pprint
import sys
import time

from PyQt4.QtCore import *

import numpy as np
import scipy as sp
from scipy.sparse import linalg
from scipy.sparse import find
from scipy.sparse import issparse
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix


from MBD_system import inverse_blockwise
from gaussian_elimination import gaussian_elimination
from MBD_system.force.force import Force


class ODE_fun(object):
    """
    classdocs
    """

    def __init__(self, MBD_system, parent=None):
        """
        DAE - Automatic generation of Ordinary Differential Equations
        :param MBD_system:
        :param parent:
        """
        #   parent
        self._parent = parent

        #   pointer to MBD system as object attribute
        self.MBD_system = MBD_system

        #   variables during integration
        self.error = False

    def evaluate_dq(self, t, q):
        """

        :param t:
        :param q:
        :return:
        """
        dy = 5. * np.exp(5. * t) * ((q - t)**2) + 1.
        return dy

    def evaluate_Jacobian(self, t, q):
        """

        :param t:
        :param q:
        :return:
        """
        dy_y = 10. * np.exp(5. * t) * (q - t)
        return dy_y

    def preprocessing(self):
        """

        :return:
        """


