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
    def __init__(self, _type, coef_of_friction_dynamic, coef_of_friction_static, parent=None):
        """
        Constructor of friction model class
        supported input for parameter _type:
        ideal
        coulomb

        """
        self._parent = parent


        self.coef_of_friction_dynamic = coef_of_friction_dynamic
        self.coef_of_friction_static = coef_of_friction_static

        #    type
        #   types:
        #   -   ideal - no friction
        #   -   friction
        self._type = _type
        if self._type == "coulomb":
            pass
        self._set_parameters_friction()


    def _set_parameters_friction(self):
        """

        :return:
        """
        self._dq_t_TOL = 1E-3


    def friction_force(self, Fn, _dq_t):
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
