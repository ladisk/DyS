"""
Created on 9. jul. 2014

@author: lskrinjar
"""
from pprint import pprint
import numpy as np


from MBD_system.fix_string import fix_string


class ContactModel(object):
    """
    classdocs
    """
    def __init__(self, _type, K=None, c_r=None, properties_dict={}, _substring="contact_model", parent=None):
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
        #   parent
        self._parent = parent

        #   sub string
        self._substring = _substring

        #   coefficient of restitution
        self.c_r = c_r

        #   initial contact velocity (normal and tangent direction)
        self._dq0_n = None
        self._dq0_t = None

        #   max penetration depth
        self._delta_max = 0.

        #   nonlinear power exponent
        self.n = 3/2.   #default is 3/2

        self.m = None
        #   length of contact (required for some types)
        self.L = None

        #   exponent for model DF
        #   list of modulus of elasticity of bodies in contact
        self.E_list = []
        self.ni_list = []

        #   status of initial contact velocity in normal and tangent direction
        self.initial_contact_velocity_calculated = False

        #   type attribute
        self._type = _type

        #   supported types
        self._types = self._supported_types()

        #   hysteresis damping model
        self._hysteresis_damping_model = None
        #   hysteresis damping models that are supported
        self._hysteresis_damping_models = ["hunt-crossley",
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

        #   if parent is specified (if exists)
        if self._parent is not None:
            #   get attributes from parent

            if hasattr(self._parent, "R0_i") and hasattr(self._parent, "R0_j"):
                self.R0_i = self._parent.R0_i
                self.R0_j = self._parent.R0_j

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

    @staticmethod
    def _supported_types():
        """

        :return:
        """
        _types = ["hertz",
                    "kelvin-voigt",
                    "hunt-crossley",
                    "herbert-mcwhannell",
                    "lee-wang",
                    "lankarani-nikravesh",
                    "gonthier et al",
                    "zhiying-qishao",
                    "flores et al",
                    "hu-guo"]
        return _types

    def _type_check(self):
        """

        :return:
        """
        if self._type.lower() in self._types:
            pass
        else:
            raise ValueError, "Contact model type not correct! Define contact model type."

    def _add_aditional_parameters(self, _dict):
        """

        :param dict:
        :return None:
        """
        for key in _dict:
            _key = key
            if hasattr(self, key):
                delattr(self, key)

            key = _key
            setattr(self, key, _dict[key])

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

        K = (4./(3.*(self.h_i + self.h_j)))*np.sqrt((self.R_i * self.R_j)/(self.R_i + self.R_j))
        return K

    def _evaluate_K_plane_sphere(self):
        """

        :return:
        """
        [self.h_i, self.h_j] = self._evaluate_h()

        K = (4/(3 * np.pi * (self.h_i + self.h_j)))*(self.R0_j**0.5)
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
        h = (1.-ni**2)/E
        return h

    def set_dq0(self, dq0_n, dq0_t):
        """
        Function saves the initial normal and tangential contact velocity at impact
        :param dq0_n    a normal contact velocity (scalar, float value)
        :param dq0_t    a tangential contact velocity (scalar, float value)
        """

        self._dq0_n = dq0_n
        self._dq0_t = dq0_t

    def reset_dq0(self):
        """

        :return:
        """
        # print "reset_dq0()"
        # self._dq0_n = None
        # self._dq0_t = None

    def contact_force(self, delta, dq_n, dq0_n=None):#dq_t, n, t
        """
        Function evaluates normal contact force and tangent friction force based on normal contact force
        :param delta:
        :param dq_n:
        :param dq_t:
        :param n:
        :param t:
        :return:
        """
        if dq0_n is not None:
            self._dq0_n = dq0_n

        #   set max penetration depth
        self._delta_max = 0.

        return self._contact_force(delta, dq_n)

    def _contact_force(self, delta, dq_n):
        """

        :param delta:
        :param dq_n:
        :return:
        """
        self._dq_n = dq_n
        if self._type == "hertz":
            Fn = self._Fn_H(delta, dq_n)
        elif self._type == "kelvin-voigt":
            Fn = self._Fn_KV(delta, dq_n)
        elif self._type == "lankarani-nikravesh":
            Fn = self._Fn_LN(delta, dq_n)
        elif self._type.lower() == "hunt-crossley":
            Fn = self._Fn_HC(delta, dq_n)
        elif self._type.lower() == "herbert-mcwhannell":
            Fn = self._Fn_HM(delta, dq_n)
        elif self._type.lower() == "lee-wang":
            Fn = self._Fn_LW(delta, dq_n)
        elif self._type.lower() == "flores et al":
            Fn = self._Fn_Flores_etal(delta, dq_n)
        elif self._type.lower() == "zhiying-qishao":
            Fn = self._Fn_ZQ(delta, dq_n)
        elif self._type.lower() == "gonthier et al":
            Fn = self._Fn_Gonthier_etal(delta, dq_n)
        elif self._type.lower() == "hu-guo":
            Fn = self._Fn_HG(delta, dq_n)
        elif self._type.lower() == "dubowsky-freudenstein":
            Fn = self._Fn_DF(delta)
        elif self._type.lower() == "esdu-78035":
            Fn = self._Fn_ESDU(delta)
        else:
            Fn = 0
            print "_type not supported: ", self._type

        self.contact = True
        if Fn < 0:
            Fn = 0.
            # self.contact = False
            # print "self._dq0_n =", self._dq0_n
            # print "self._dq_n =", self._dq_n
            # print self._Fn_LN(delta, dq_n)
            # print "Fn is set to 0, otherwise it would be negative!"
            # print "step =", self._parent._step

        return Fn

    def _Fn_H(self, delta, _dq_n):
        """
        Hertz contact force
        :param delta:
        :param _dq_n:
        :return:
        """
        _Fn = self.K * abs(delta)**self.n
        return _Fn

    def _Fn_KV(self, delta, _dq_n):
        """
        Kelvin-Voigt contact force
        :return:
        """
        #   compression phase
        if _dq_n >= 0:
            _Fn = self.K * abs(delta)
            return _Fn
        #   restitution phase
        elif _dq_n < 0:
            _Fn = self.c_r * self.K * abs(delta)
            return _Fn

    def _Fn(self, delta, _dq_n):
        """
        Hertz contact force model
        :return:
        """
        F = self.K * abs(delta)**self.n * (1 + self.ksi * (_dq_n / self._dq0_n))#(self.ksi * (abs(delta)**self.n) * (_dq_n / self._dq0_n) * self.K)
        return F

    def _Fn_LN(self, delta, _dq_n):
        """
        Lankarani-Nikravesh contact force model
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self._hysteresis_damping_factor_LN()
        #    normal force
        _Fn = self._Fn(delta, _dq_n)
        return _Fn

    def _Fn_HC(self, delta, _dq_n):
        """
        Hunt and Crossley contact force model
        :param delta:
        :param _dq_n:
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self._hysteresis_damping_factor_HC()
        #    normal force
        _Fn = self._Fn(delta, _dq_n)
        return _Fn

    def _Fn_HM(self, delta, _dq_n):
        """
        Herbert-McWhannell contact force model
        :param delta:
        :param _dq_n:
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self._hysteresis_damping_factor_HM()
        #    normal force
        _Fn = self._Fn(delta, _dq_n)
        return _Fn

    def _Fn_LW(self, delta, _dq_n):
        """
        Lee-Wang contact force model
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self._hysteresis_damping_factor_LW()
        #   contact force
        _Fn = self._Fn(delta, _dq_n)
        return _Fn

    def _Fn_Flores_etal(self, delta, _dq_n):
        """
        Flores et al. contact force model
        :param delta:
        :param _dq_n:
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self._hysteresis_damping_factor_Flores_etal()
        #   contact force
        _Fn = self._Fn(delta, _dq_n)
        return _Fn

    def _Fn_ZQ(self, delta, _dq_n):
        """
        Zhiying-Qishao contact force model
        :param delta:
        :param _dq_n:
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self._hysteresis_damping_factor_ZQ()
        #   contact force
        _Fn = self._Fn(delta, _dq_n)
        return _Fn

    def _Fn_Gonthier_etal(self, delta, _dq_n):
        """
        Gonthier et al. contact force model
        :param delta:
        :param _dq_n:
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self._hysteresis_damping_factor_Gonthier_etal()
        #   contact force
        _Fn = self._Fn(delta, _dq_n)
        return _Fn

    def _Fn_HG(self, delta, _dq_n):
        """
        Hu-Guo contact force model
        :param delta:
        :param _dq_n:
        :return:
        """
        #   hysteresis damping factor
        self.ksi = self._hysteresis_damping_factor_HG()

        #   contact force
        _Fn = self._Fn(delta, _dq_n)
        return _Fn

    def hysteresis_damping_factor(self):
        """
        Function evaluates hysteresis damping factor based on type of contact model used only for testing not during simulation
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
        return ksi * (self.K/self._dq0_n)

    def _hysteresis_damping_factor_LN(self):
        """
        Lankarani-Nikravesh
        :return:
        """
        ksi = (3./4.)*(1. - self.c_r ** 2)
        return ksi

    def _hysteresis_damping_factor_HC(self):
        """
        Hunt-Crossley
        :return:
        """
        ksi = (3/2.)*(1-self.c_r)
        return ksi

    def _hysteresis_damping_factor_HM(self):
        """

        :return:
        """
        ksi = (6.*(1. - self.c_r))/((2.*self.c_r-1.)**2 + 3.)
        return ksi

    def _hysteresis_damping_factor_LW(self):
        """

        :return:
        """
        ksi = (3./4.)*(1. - self.c_r)
        return ksi

    def _hysteresis_damping_factor_Flores_etal(self):
        """
        Flores et al
        :return:
        """
        ksi = (8*(1-self.c_r))/(5*self.c_r)
        return ksi

    def _hysteresis_damping_factor_ZQ(self):
        """

        :return:
        """
        ksi = ((3/4.)*(1-self.c_r**2)*np.exp(2*(1-self.c_r)))
        return ksi

    def _hysteresis_damping_factor_Gonthier_etal(self):
        """
        Gonthier et al
        :return:
        """
        #   dimensionless factor approximated by
        d = 1-self.c_r**2

        #   hysteresis damping factor
        ksi = (d/self.c_r)
        return ksi

    def _hysteresis_damping_factor_HG(self):
        """

        :return:
        """
        ksi = ((3/2.)*(1-self.c_r)/(2*self.c_r))
        return ksi


if __name__ == "__main__":
    models = []
    K = 19224012713.2
    cr = .7
    delta = -2E-6
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
    model_type = "hertz"
    # model_type = "dubowsky-freudenstein"
    # model_type = "esdu-78035"
    # model_type = "lankarani-nikravesh"
    #   material properties for aluminium
    E = 69E+9
    ni = 0.35

    _model = ContactModel(_type=model_type, K=K, c_r=cr)
    _model.m = 3
    _model.L = 4.E-3
    _model.h_i = _model._h(E, ni)
    # print "_model.h_i =", _model.h_i
    _model.h_j = _model._h(E, ni)
    _model.R_i = 0.011
    #   contact for flat section
    _model.R_i = 1E+2
    _model.R_j = 0.01
    _model._maxiter = 50
    _model._tol = 1.e-06
    # _model.K = 1E+6
    _model._dq0_n = -0.0284053065794
    delta = -1.90798765398e-06
    # pprint(vars(_model))
    #   contact force
    # pprint(vars(_model))
    Fn = _model.contact_force(delta, dq_n=-0.0284053065794)
    print "Fn(sol) =", Fn