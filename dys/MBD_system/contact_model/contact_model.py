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


class ContactModel(object):
    '''
    classdocs
    '''
    def __init__(self, _type, coef_of_restitution, parent=None):
        """
        Constructor of contact model class
        type:
        """
        self._parent = parent

        self.coef_of_restitution = coef_of_restitution

        self.initial_contact_velocity_calculated = False
        #   type
        #   types:
        #   -   kelvin-voigt
        #   -   lankarani-nikravesh
        self._type = _type

        if self._type == "lankarani-nikravesh":
            self.__set_parameters_LN()
        elif self._type == "kelvin-voigt":
            self.__set_parameters_LN()
        elif self._type == "hertz":
            self.__set_parameters_H()
        else:
            raise ValueError, "Contact model type not defined! Define contact model type."


    def _h(self, E, ni):
        """

        :return:
        """
        h = (1-ni**2)/(np.pi * E)
        return h


    def __set_parameters_H(self):
        """

        :return:
        """
        #   module of elasticity
        E_i = 2.1E+5
        E_j = 2.1E+5
        #   ni module
        ni_i = 0.3
        ni_j = 0.3
        #   value of h
        h_i = self._h(E_i, ni_i)
        h_j = self._h(E_j, ni_j)


        self.K = (4/(3. * np.pi * (h_i+h_j))) * ((self._parent.R_i*self._parent.R_j) / (abs(self._parent.R_i-self._parent.R_j)))**(0.5)


    def __set_parameters_LN(self):
        """

        :return:
        """
        #   generalized stiffness
        self.K = 6.0E+14
        #   exponent
        self.n = 1.5
        #   coefficient of restitution
        self.coef_of_restitution = 0.5


    def __set_parameters_KV(self):
        """

        :return:
        """
        #   generalized stiffness
        self.K = 5E+12


    def set_dq0(self, dq0_n, dq0_t):
        """
        
        """
        self._dq0_n = dq0_n
        self._dq0_t = dq0_t


    def contact_force(self, delta, dq_n, dq_t, n, t):
        """
        Function calculates normal contact force and tangent friction force based on normal contact force
        :param delta:
        :param dq_n:
        :param dq_t:
        :param n:
        :param t:
        :return:
        """
        if self._type == "lankarani-nikravesh":
            Fn = self.__contact_force_LN(delta, dq_n)
        elif self._type == "kelvin-voigt":
            Fn = self.__contact_force_KV(delta, dq_n)
        elif self._type == "hertz":
            Fn = self.__contact_force_H(delta, dq_n)
        else:
            raise ValueError, "Contact model type not defined! Define contact model type."


        return Fn


    def __contact_force_H(self, delta, _dq_n):
        """

        :param delta:
        :param _dq_n:
        :return:
        """
        _Fn = self.K * (abs(delta))
        return _Fn


    def __contact_force_KV(self, delta, _dq_n):
        """

        :return:
        """
        #    normal force
        if _dq_n > 0:
            _Fn = self.K * (abs(delta))
            return _Fn
        elif _dq_n < 0:
            _Fn = self.coef_of_restitution * self.K * (abs(delta))
            return _Fn


    def __contact_force_LN(self, delta, _dq_n):
        """

        :return:
        """
        #    normal force
        _Fn = self.K * (abs(delta) ** self.n) * (1+0.75*(1 - self.coef_of_restitution ** 2)*(_dq_n/self._dq0_n))
        return _Fn