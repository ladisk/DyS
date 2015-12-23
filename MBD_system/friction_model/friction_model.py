'''
Created on 9. jul. 2014

@author: lskrinjar
'''

import numpy as np


try:
    from ..MBD_system import *
    from ..cad2cm_lcs import cad2cm_lcs
    from ..force_matrix import Force_Q_e_matrix
except:
    None


class FrictionModel(object):
    '''
    classdocs
    '''
    def __init__(self, _type, coef_of_friction_kinematic=None, coef_of_friction_dynamic=None, coef_of_friction_static=None, parent=None):
        """
        Constructor of friction model class
        supported input for parameter
        _type:
            Ideal (no friction)
            Coulomb (friction)
            Viscous (friction)
            Stribeck (friction)
            Ambrosio (friction)
            LuGre (friction)
            Dahl (friction)
            Bristle (friction)
            Reset integrator (friction)
            Karnopp (friction)
            Bliman-Sorine (friction)
            Leuven (friction)

        """
        self._parent = parent

        #   type
        self._type = _type

        #   coefficient of static friction
        self.coef_of_friction_static = coef_of_friction_static

        #   coulomb friction model
        #   coefficient of kinematic friction
        self.coef_of_friction_kinematic = coef_of_friction_kinematic

        #   viscous friction model
        #   coefficient of viscous friction
        self.sigma_v = None

        #   stiction friction model
        #   sliding speed coefficient
        self.v_s = None
        #   gradient of friction decay
        self.gamma = None

        #   ambrosio friction model
        self.c_d = None
        self.v0 = None
        self.v1 = None

        if self._type == "coulomb":
            pass
        self._set_parameters_friction()

    def _set_parameters_friction(self):
        """

        :return:
        """
        self._dq_t_TOL = 1E-3

    def friction_force_testing(self, Fn, _dq_t):
        """

        :param Fn:
        :param _dq_t:
        :return:
        :param Ft: tangent force scalar value
        """
        #   tangent force
        if abs(_dq_t) < self._dq_t_TOL:
            coef_of_friction = self.coef_of_friction_static
        else:
            coef_of_friction = self.coef_of_friction_dynamic
            

        # print "coef_of_friction =", coef_of_friction
        Ft = Fn * coef_of_friction
        
        return Ft

    def friction_force(self, Fn, _dq_t):
        """

        :return:
        """

        if self._type.lower() == "ideal":
            Ft = 0
        elif self._type.lower() == "coulomb":
            Ft = self.__friction_force_coulomb(Fn, _dq_t)
        elif self._type.lower() == "coulomb":
            Ft = self.__friction_force_coulomb(Fn, _dq_t)
        elif self._type.lower() == "viscous":
            Ft = self.__friction_force_viscous(Fn, _dq_t)
        elif self._type.lower() == "stribeck":
            Ft = self.__friction_force_stribeck(Fn, _dq_t)
        elif self._type.lower() == "ambrosio":
            Ft = self.__friction_force_ambrosio(Fn, _dq_t)
        elif self._type.lower() == "lugre":
            Ft = self.__friction_force_lugre(Fn, _dq_t)
        elif self._type.lower() == "coulomb":
            Ft = self.__friction_force_coulomb(Fn, _dq_t)
        elif self._type.lower() == "coulomb":
            Ft = self.__friction_force_coulomb(Fn, _dq_t)
        elif self._type.lower() == "coulomb":
            Ft = self.__friction_force_coulomb(Fn, _dq_t)
        elif self._type.lower() == "coulomb":
            Ft = self.__friction_force_coulomb(Fn, _dq_t)
        else:
            raise ValueError, "Contact model type not defined! Define contact model type."
        return Ft

    def __friction_force_coulomb(self, Fn, _dq_t):
        """
        Coulomb friction model
        :param Fn:
        :param _dq_t:
        :return:
        """
        if _dq_t == 0:
            _Ft = 0#Q_e
        else:
            _Ft = -self.coef_of_friction_kinematic * Fn * np.sign(_dq_t)
        return _Ft

    def __friction_force_viscous(self, Fn, _dq_t):
        """
        Viscous friction model with coulomb friction
        :param Fn:
        :param _dq_t:
        :return:
        """
        Fc = self.__friction_force_coulomb(Fn, _dq_t)
        _Ft = Fc + np.sign(Fc) * self.sigma_v * _dq_t
        return _Ft

    def __friction_force_stribeck(self, Fn, _dq_t):
        """
        Stribeck friction model
        :param Fn:
        :param _dq_t:
        :return:
        """
        Fc = self.__friction_force_coulomb(Fn, _dq_t)
        Fs = -self.coef_of_friction_static * Fn * np.sign(_dq_t)
        _Ft = (Fc+(Fs-Fc)*np.exp((-abs(_dq_t)/self.v_s)**self.i))*np.sign(_dq_t) + self.sigma_v * _dq_t
        return _Ft

    def __friction_force_ambrosio(self, Fn, _dq_t):
        """
        Ambrosio friction model
        :param Fn:
        :param _dq_t:
        :return:
        """
        self.c_d = self._evaluate_c_d(_dq_t)

        _Ft = -self.coef_of_friction_kinematic * self.c_d * Fn * np.sign(_dq_t)
        return _Ft

    def _evaluate_c_d(self, _dq_t):
        """

        :param Fn:
        :param _dq_t:
        :return:
        """
        if _dq_t <= self.v0:
            c_d = 0

        if self.v0 <= _dq_t <= self.v1:
            c_d = (_dq_t - self.v0)/(self.v1 - self.v0)

        if _dq_t >= self.v1:
            c_d = 1

        return c_d

    def __friction_force_lugre(self, Fn):
        """
        LuGre friction model
        :return:
        """


        _Ft = self.coef_of_friction * Fn
        return _Ft

    def _evaluate_lugre_coef_of_friction(self):
        """

        :return:
        """
