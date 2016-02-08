"""

created by: lskrinjar
date of creation: 01/02/2016
time of creation: 10:11
"""
import itertools
from OpenGL.GL import *
import numpy as np
from pprint import pprint
import sympy as sp
from sympy import Symbol
from sympy.parsing.sympy_parser import parse_expr


from MBD_system.MBD_system_items import MotionItem
from MBD_system.q2R_i import q2R_i
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2theta_i import q2theta_i


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

        #   positions and rotations
        self.q = [self.Rx, self.Ry, self.theta]
        #   velocities (translational and rotational)
        self.dq = [self.dRx, self.dRy, self.dtheta]

        #   generalized displacements and velocities vector of a body
        self._q = [self.Rx, self.Ry, self.theta, self.dRx, self.dRy, self.dtheta]

        self._q_eval = []

        #   attribute names
        self.__q_names = ["Rx", "Ry", "theta"]
        self.__dq_names = ["dRx", "dRy", "dtheta"]

        #   get q vector of current MBD system at joint object initialization
        self.q0 = self._parent._parent.get_q()
        
        #    use symbolic input (as string) to define functions of time
        self.t = Symbol('time')

        #   evaluate number of constraints
        self.__number_of_constraints()

    def __number_of_constraints(self):
        """

        :return:
        """
        for q_i, q_name in zip(self.q, self.__q_names):
            if q_i is not None:
                self._q_eval.append(q_name)

        self.C_size = len(self._q_eval)

    def evaluate_C(self, q, t=0):
        """
        Function evaluates constraint equations on positions level
        :param q:   vector of MBD system
        :param q:   time
        """
        if q is None:
            q = self._parent._parent.get_q()

        if t is False:
            t = 0

        #   predefine vector
        C = np.zeros(self.C_size)
        # print "C(initial) =", C
        j = 0
        # print "self.__q_names =", self.__q_names
        for i, (q_i, q_name) in enumerate(zip(self.q, self.__q_names)):
            # print "q_i =", q_i
            if q_name == "Rx":
                dC = (q2R_i(q, self.body_id) - q2R_i(self.q0, self.body_id))[0]
            if q_name == "Ry":
                dC = (q2R_i(q, self.body_id) - q2R_i(self.q0, self.body_id))[1]
            if  q_name == "theta":
                dC = q2theta_i(q, self.body_id) - q2theta_i(self.q0, self.body_id)

            if q_i is not None:
                _str = getattr(self, q_name)
                #   evaluate expression
                val = eval(_str, {}, {"time":t})
                _C = dC - val
                # print "_C =", _C
                C[j] = _C
                j += 1
        # print "C(final) =", C
        return C

    def evaluate_C_t(self, q=None, t=None):
        """
        Function evaluates a time derivative of constraint equation
        """
        C_t = np.zeros(self.C_size)
        j = 0
        for i, (q_i, q_name) in enumerate(zip(self.q, self.__q_names)):
            if q_i is not None:
                _exp = parse_expr(getattr(self, q_name))
                #   evaluate derivative
                _dexpr = sp.diff(_exp, self.t)
                C_t[j] = _dexpr
                j += 1

        return C_t

    def evaluate_C_q(self, q):
        """

        :param q:
        :return:
        """
        C_q = np.zeros(len(self.q0)/2)
        j = 0
        for i, (q_i, q_name) in enumerate(zip(self.q, self.__q_names)):
            if q_i is not None:
                C_q[3*self.body_id + i] = 1
                j += 1

        return C_q

    def evaluate_Q_d(self, q, t):
        """

        :param q:
        :param t:
        :return:
        """
        return np.zeros(self.C_size)
