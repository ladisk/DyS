"""

created by: lskrinjar
date of creation: 01/02/2016
time of creation: 10:11
"""
import itertools
from OpenGL.GL import *
import numpy as np
from sympy import Symbol
from sympy.parsing.sympy_parser import parse_expr

from pprint import pprint

from MBD_system.MBD_system_items import MotionItem


class Motion(MotionItem):
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, name, body_id, Rx=None, Ry=None, theta=None, dRx=None, dRy=None, dtheta=None, parent=None):
        """
        Constructor of a class Motion that can be used in kinematic analysis of the system

        """
        super(Motion, self).__init__(name, parent)

        #   motion id
        self._id = self.__id.next()

        #   parent
        self._parent = parent

        #   name
        self._name = name

        #   body id - pointer
        self.body_id = body_id

        #   generalized coordinates of a body
        #   explicit function of time
        #   Rx
        self.Rx = Rx
        #   Ry
        self.Ry = Ry
        #   theta
        self.theta = theta
        #   dRx
        self.dRx = dRx
        #   dRy
        self.dRy = dRy
        #   dtheta
        self.dtheta = dtheta

        #   generalized displacements and velocites vector of a body
        self.q = [self.Rx, self.Ry, self.theta, self.dRx, self.dRy, self.dtheta]

    def _evaluate(self, t, q):
        """
        Function evaluates coordinate at time
        :param t:   time
        :param q:   vector of absolute coordinates and velocities of the system
        :return:
        """
        for q_i in self.q:
            if q_i is not None:
                pass
