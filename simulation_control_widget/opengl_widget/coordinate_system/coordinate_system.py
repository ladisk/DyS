"""
Created on 9. jul. 2015

@author: reti-luka
"""
import sys
import os
from pprint import pprint
import itertools
import logging
import numpy as np
from PyQt4 import QtCore, QtGui, QtOpenGL


from simulation_control_widget.opengl_widget.marker.marker import Marker


class CoordinateSystem(Marker):
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, parent=None):
        """

        :param parent:
        :return:
        """
        #    number
        self._parent = parent

        #   coordinate system id
        self.coordinate_system_id = self.__id.next()

        if self._parent is None:
            self._name = "GCS"
        else:
            self._name = "CS_"+str(self._parent.body_id)

        super(CoordinateSystem, self).__init__(node=np.array([0., 0., 0.], dtype=np.float32), visible=True, scale=1E-2, parent=parent)