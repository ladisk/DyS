"""
Created on 9. jul. 2014

@author: lskrinjar
"""
import time
import numpy as np
from matplotlib import pyplot as plt
from pprint import pprint
import scipy
import scipy.optimize


from MBD_system.fix_string import fix_string

class ContactModel(object):
    """
    classdocs
    """
    def __init__(self, _type, K=None, c_r=None, properties_dict=[], parent=None):
        """
        Constructor of contact model class
        Contact force models are built based on document:
        Compliand contact force models in multibody dynamics: Evolution of the Hertz theory, Machado et al, 2012
        :param _type:   supported types for contact model
                        1 Hertz
                        2 Kelvin-Voigt
                        3 Hunt-Crossley
                        4 Herbert-McWhannell
                        5 Lee-Wang
                        6 Lankarani-Nikravesh
                        7 Gonthier et al.
                        8 Zhiying-Qishao
                        9 Flores et al.
                        10 Hu-Guo
                        11 Dubowsky-Freudenstein
                        12 ESDU-78035 (todo)
        :param  K       generalized contact stiffness, can be defined by user input
        :param  c_r     coefficient of restitution
        :parem  parent  pointer to the parent object
        """
        self._parent = parent

        #   sub string
        self._substring = "contact_model"

        #   if parent is specified (if exists)
        if self._parent is not None:
            #   get attributes from parent
            if hasattr(self._parent, "R_i") and hasattr(self._parent, "R_j"):
                self.R_i = self._parent.R_i
                self.R_j = self._parent.R_j

            if hasattr(self._parent, "K"):
                self.K = self._parent.K

            if self._parent._parent is not None and K is None:
                if self._parent._type.lower() == "contact sphere-sphere":#contact sphere-sphere
                    self.K = self._evaluate_K_sphere_sphere()
                if self._parent._type.lower() == "contact plane-sphere":
                    self.K = self._evaluate_K_plane_sphere()

        #   if parent is None
        else:
            self.K = K
            self.R_i = None
            self.R_j = None

        #   coefficient of restitution
        self.c_r = c_r

        #   nonlinear power
        self.n = 3/2.   #default is 3/2

        #   length of contact (required for some types)
        self.L = None

        #   exponent for model DF
        self.m = None

        #   list of modulus of elasticity of bodies in contact
        self.E_list = []
        self.ni_list = []

        #   contact force value from previous evaluation
        self._Fn_last = 1E-3

        #   number of iterations for implicit function solver
        self._maxiter = 100
        self._tol = 1e-06

        #   status of initial contact velocity in normal and tangent direction
        self.initial_contact_velocity_calculated = False

        #   type attribute
        self._type = _type
        self._types = ["hertz",
                        "kelvin-voigt",
                        "hunt-crossley",
                        "herbert-mcwhannell",
                        "lee-wang",
                        "lankarani-nikravesh",
                        "gonthier et al",
                        "zhiying-qishao",
                        "flores et al",
                        "hu-guo"]

        #   set parameters by type
        if isinstance(self._type, basestring) or self._type is None:
            if self._type is None:
                self._type = "hertz"

        #   set additional properties to object from dictionary
        self.properties = properties_dict
        if self.properties is not []:
            self._add_aditional_parameters(self.properties)

    def _type_check(self):
        """

        :return:
        """
        print
        if self._type.lower() in self._types:
            pass
        else:
            raise ValueError, "Contact model type not correct! Define contact model type."

    def _add_aditional_parameters(self, dict):
        """

        :param dict:
        :return None:
        """
        for key in dict:
            if hasattr(self, key):
                delattr(self, key)

            setattr(self, key, dict[key])

    def setType(self, type):
        """

        :param type:
        :return:
        """
        self._type = fix_string(type.title())

    def _evaluate_K_sphere_sphere(self):
        """

        :return:
        """
        [self.h_i, self.h_j] = self._evaluate_h()

        K = (4/(3*(self.h_i + self.h_j)))*np.sqrt((self._parent.R_i * self._parent.R_j)/(self._parent.R_i + self._parent.R_j))
        return K

    def _evaluate_K_plane_sphere(self):
        """

        :return:
        """
        [self.h_i, self.h_j] = self._evaluate_h()

        K = (4/(3 * np.pi * (self.h_i + self.h_j)))*(self._parent.R0_j**0.5)
        return K

    def _evaluate_h(self):
        """

        :return:
        """
        self.h_list = []
        for body_id in self._parent.body_id_list:
            body = self._parent._parent._parent.bodies[body_id]

            _h = self._h(body.module_of_elasticity, body.poisson_ratio)

            self.h_list.append(_h)

        [h_i, h_j] = self.h_list

        return h_i, h_j

    def _h(self, E, ni):
        """

        :return:
        """
        h = (1-ni**2)/E
        return h

    def set_dq0(self, dq0_n, dq0_t):
        """
        Function saves the initial normal and tangential contact velocity at impact
        :param dq0_n    a normal contact velocity (scalar, float value)
        :param dq0_t    a tangential contact velocity (scalar, float value)
        """
        self._dq0_n = dq0_n
        self._dq0_t = dq0_t

    def contact_force(self, delta, dq_n):#dq_t, n, t
        """
        Function calculates normal contact force and tangent friction force based on normal contact force
        :param delta:
        :param dq_n:
        :param dq_t:
        :param n:
        :param t:
        :return:
        """
        if self._type == "hertz":
            Fn = self.__Fn_H(delta, dq_n)
        elif self._type == "kelvin-voigt":
            Fn = self.__Fn_KV(delta, dq_n)
        elif self._type == "lankarani-nikravesh":
            Fn = self.__Fn_LN(delta, dq_n)
        elif self._type.lower() == "hunt-crossley":
            Fn = self.__Fn_HC(delta, dq_n)
        elif self._type.lower() == "herbert-mcwhannell":
            Fn = self.__Fn_HM(delta, dq_n)
        elif self._type.lower() == "lee-wang":
            Fn = self.__Fn_LW(delta, dq_n)
        elif self._type.lower() == "flores et al":
            Fn = self.__Fn_Flores_etal(delta, dq_n)
        elif self._type.lower() == "zhiying-qishao":
            Fn = self.__Fn_ZQ(delta, dq_n)
        elif self._type.lower() == "gonthier et al":
            Fn = self.__Fn_Gonthier_etal(delta, dq_n)
        elif self._type.lower() == "hu-guo":
            Fn = self.__Fn_HG(delta, dq_n)
        elif self._type.lower() == "dubowsky-freudenstein":
            Fn = self.__Fn_DF(delta)
        elif self._type.lower() == "esdu-78035":
            Fn = self.__Fn_ESDU(delta)
        else:
            raise ValueError, "Contact model type not defined! Define contact model type."
        # print "Fn =", Fn, "delta =", delta, "dqn =", dq_n#, "type =", self._type
        return Fn

    def __Fn_H(self, delta, _dq_n):
        """
        Hertz contact force
        :param delta:
        :param _dq_n:
        :return:
        """
        _Fn = -self.K * abs(delta)**self.n
        return _Fn

    def __Fn_KV(self, delta, _dq_n):
        """
        Kelvin-Voigt contact force
        :return:
        """
        #   compression phase
        if _dq_n < 0:
            _Fn = -self.K * abs(delta)
            return _Fn
        #   restitution phase
        elif _dq_n >= 0:
            _Fn = -self.c_r * self.K * abs(delta)
            return _Fn

    def __Fn(self, delta, _dq_n):
        """
        Hertz contacr force model
        :return:
        """
        F = -self.K * abs(delta)**self.n - self.ksi * (abs(delta)**self.n) * _dq_n
        return F

    def __Fn_LN(self, delta, _dq_n):
        """
        Lankarani-Nikravesh contact force model
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self._hysteresis_damping_factor_LN()
        #    normal force
        _Fn = self.__Fn(delta, _dq_n)
        return _Fn

    def __Fn_HC(self, delta, _dq_n):
        """
        Hunt and Crossley contact force model
        :param delta:
        :param _dq_n:
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self._hysteresis_damping_factor_HC()
        #    normal force
        _Fn = self.__Fn(delta, _dq_n)
        return _Fn

    def __Fn_HM(self, delta, _dq_n):
        """
        Herbert-McWhannell contact force model
        :param delta:
        :param _dq_n:
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self._hysteresis_damping_factor_HM()
        #    normal force
        _Fn = self.__Fn(delta, _dq_n)
        return _Fn

    def __Fn_LW(self, delta, _dq_n):
        """
        Lee-Wang contact force model
        :return:
        """
        #   hysteresis damping factor
        self.ksi = (3/4.)*(1-self.c_r) * (self.K/self._dq0_n)
        #   contact force
        _Fn = self.__Fn(delta, _dq_n)
        return _Fn

    def __Fn_Flores_etal(self, delta, _dq_n):
        """
        Flores et al. contact force model
        :param delta:
        :param _dq_n:
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self._hysteresis_damping_factor_Flores_etal()
        #   contact force
        _Fn = self.__Fn(delta, _dq_n)
        return _Fn

    def __Fn_ZQ(self, delta, _dq_n):
        """
        Zhiying-Qishao contact force model
        :param delta:
        :param _dq_n:
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self._hysteresis_damping_factor_ZQ()
        #   contact force
        _Fn = self.__Fn(delta, _dq_n)
        return _Fn

    def __Fn_Gonthier_etal(self, delta, _dq_n):
        """
        Gonthier et al. contact force model
        :param delta:
        :param _dq_n:
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self.__Fn_Gonthier_etal()
        #   contact force
        _Fn = self.__Fn(delta, _dq_n)
        return _Fn

    def __Fn_HG(self, delta, _dq_n):
        """
        Hu-Guo contact force model
        :param delta:
        :param _dq_n:
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self._hysteresis_damping_factor_HG()
        _Fn = self.__Fn(delta, _dq_n)
        return _Fn

    def hysteresis_damping_factor(self):
        """
        Function evaluates hysteresis damping factor based on type of contact model
        :return:
        """
        if self._type.lower() == "hertz":
            pass
        elif self._type.lower() == "kelvin-voigt":
            pass
        elif self._type.lower() == "lankarani-nikravesh":
            ksi = self._hysteresis_damping_factor_LN()
        elif self._type.lower() == "hunt-crossley":
            ksi = self._hysteresis_damping_factor_HC()
        elif self._type.lower() == "herbert-mcwhannell":
            ksi = self._hysteresis_damping_factor_HM()
        elif self._type.lower() == "lee-wang":
            ksi = self._hysteresis_damping_factor_LW()
        elif self._type.lower() == "flores et al":
            ksi = self._hysteresis_damping_factor_Flores_etal()
        elif self._type.lower() == "zhiying-qishao":
            ksi = self._hysteresis_damping_factor_ZQ()
        elif self._type.lower() == "gonthier et al":
            ksi = self._hysteresis_damping_factor_Gonthier_etal()
        elif self._type.lower() == "hu-guo":
            ksi = self._hysteresis_damping_factor_HG()
        else:
            raise ValueError, "Contact model type not defined! Define contact model type."
        return ksi

    def _hysteresis_damping_factor_LN(self):
        """
        Lankarani-Nikravesh
        :return:
        """
        ksi = (3/4.)*(1. - self.c_r ** 2) * (self.K/self._dq0_n)
        return ksi

    def _hysteresis_damping_factor_HC(self):
        """
        Hunt-Crossley
        :return:
        """
        ksi = (3/2.)*(1-self.c_r) * (self.K/self._dq0_n)
        return ksi

    def _hysteresis_damping_factor_HM(self):
        """

        :return:
        """
        ksi = (6*(1-self.c_r))/((2*self.c_r-1)**2 + 3) * (self.K/self._dq0_n)
        return ksi

    def _hysteresis_damping_factor_LW(self):
        """

        :return:
        """
        ksi = (3/4.)*(1-self.c_r) * (self.K/self._dq0_n) * (self.K/self._dq0_n)
        return ksi

    def _hysteresis_damping_factor_Flores_etal(self):
        """
        Flores et al
        :return:
        """
        ksi = (8*(1-self.c_r))/(5*self.c_r) * (self.K/self._dq0_n)
        return ksi

    def _hysteresis_damping_factor_ZQ(self):
        """

        :return:
        """
        ksi = ((3/4.)*(1-self.c_r**2)*np.exp(2*(1-self.c_r))) * (self.K/self._dq0_n)
        return ksi

    def _hysteresis_damping_factor_Gonthier_etal(self):
        """
        Gonthier et al
        :return:
        """
        #   dimensionless factor approximated by
        d = 1-self.c_r**2

        #   hysteresis damping factor
        ksi = (d/self.c_r) * (self.K/self._dq0_n)
        return ksi

    def _hysteresis_damping_factor_HG(self):
        """

        :return:
        """
        ksi = ((3/2.)*(1-self.c_r)/(2*self.c_r)) * (self.K/self._dq0_n)
        return ksi

if __name__ == "__main__":
    models = []
    K = 19224012713.2
    cr = .7
    delta = -2E-5
    dq_n = -.01
    dq0 = -1

    # model_types = ["hertz", "kelvin-voigt", "lankarani-nikravesh", "hunt-crossley", "herbert-mcwhannell", "lee-wang", "flores et al", "zhiying-qishao", "gonthier et al"]
    #
    # for i, model_type in enumerate(model_types):
    #     _model = ContactModel(_type=model_type, K=K, c_r=cr)
    #     _model.set_dq0(dq0, 0)
    #
    #     #   contact force
    #     Fn = _model.contact_force(delta, dq_n)
    #
    #     models.append(_model)
    #     if True: #i == 0 or i == 2 or i == 3 or i == 4 or i == 5 or i == 6:
    #         print "------------------------------"
    #         print i,
    #         print "type =", model_type,
    #         print "Fn =", Fn,
    #         if hasattr(_model, "ksi"):
    #             print "ksi =", _model.ksi
    #         else:
    #             print "ksi =", np.NaN

    #   testing implicit contact function solver
    model_type = "dubowsky-freudenstein"
    model_type = "esdu-78035"
    E = 2.1E+9
    ni = 0.3
    _model = ContactModel(_type=model_type, K=K, c_r=cr)
    _model.m = 3
    _model.L = 4.E-3
    _model.h_i = _model._h(E, ni)
    # print "_model.h_i =", _model.h_i
    _model.h_j = _model._h(E, ni)
    _model.R_i = 0.011
    _model.R_j = 0.01
    _model._maxiter = 50
    _model._tol = 1.e-06

    _model._Fn0 = 1.e-6
    # pprint(vars(_model))
    #   contact force
    Fn = _model.contact_force(delta, dq_n=0)
    print "Fn(sol) =", Fn