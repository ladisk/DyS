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

    def __init__(self, name, body_id, q0=None, q1=None, q2=None, q3=None, q4=None, q5=None, parent=None):
        """
        Constructor
        """
        super(Motion, self).__init__(name, parent)

        #   parent
        self._parent = parent

        #   name
        self._name = name

        #   body id - pointer
        self.body_id = body_id

        #   generalized coordinates of a body
        #   explicit function of time
        #   R_x
        self.q0 = q0
        #   R_y
        self.q1 = q1
        #   theta
        self.q2 = q2
        #   dR_x
        self.q3 = q3
        #   dR_y
        self.q4 = q4
        #   dtheta
        self.q5 = q5

        #   generalized displacements and velocites vector of a body
        self.q = [self.q0, self.q1, self.q2, self.q3, self.q4, self.q5]

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
