"""
Created on 1. jan. 2016

@author: lskrinjar
"""
from pprint import pprint
import copy
import matplotlib as mpl
import numpy as np
import scipy
import scipy.optimize
from matplotlib import pyplot as plt
import xlsxwriter


from MBD_system.contact_model.contact_model import ContactModel


class ContactModelCylinder(ContactModel):
    """
    classdocs
    """


    def __init__(self, _type, K=0., c_r=0., properties_dict={}, _substring="contact_model", parent=None):
        """
        If internal type of cylindrical contact is used:
        Body i - hole
        Body j - pin

        :param _type:   supported types for cylindrical contact model
                        1 Hertz
                        2 Johnson
                        3 Radzimovsky
                        4 Goldsmith
                        5 Dubowsky-Freudenstein
                        6 ESDU-78035
                        7 Liu etal
                        8 Lankarani-Nikravesh
                        9 Pereira etal
                       10 Persson
        :return:
        """
        super(ContactModelCylinder, self).__init__(_type, K=None, c_r=None, properties_dict=properties_dict, parent=parent)

        #   contact location, options:
        #   internal (default)
        #   external
        self.contact_location = "internal"
        
        #   sub string
        self._substring = _substring

        #   equivalent contact modulus of elasticity
        self._evaluate_E()

        #   reduced contact radius
        self._evaluate_dR()

        #   relative curvature of the contact (only used for Radzimovsky type)
        self._evaluate_R()

        #   length of contact
        if self._parent is None:
            self.L = None
        else:
            self.L = self._parent.L

        #   exponential value
        if hasattr(self._parent, "m"):
            self.m = self._parent.m
        else:
            self.m = None

        #   exponent
        self.n = 1.5

        #   iteration settings
        self._maxiter = 2000
        self._tol = 1E-2
        self._Fn_last = 1E-2

        #   type attribute
        #   persson not implemented
        self._type = _type
        self.__types = copy.copy(self._types)
        self._types = ["hertz cylinder",                     #done
                    "johnson",                  #done
                    "radzimovsky",              #done
                    "goldsmith",                #done
                    "dubowsky-freudenstein",    #done
                    "esdu-78035",               #done
                    "liu et al",                #done
                    "lankarani-nikravesh cylinder",      #done
                    "pereira et al"]            #done

        #   set parameters by type
        if isinstance(self._type, basestring) or self._type is None:
            if self._type is None:
                self._type = "hertz"
            else:
                self._type_check()

        #   set additional properties to object from dictionary
        self.properties = properties_dict
        if self.properties is not []:
            self._add_aditional_parameters(self.properties)
            
    @staticmethod
    def _supported_types():
        """

        :return:
        """
        _types = ["hertz cylinder",                     #done
                    "johnson",                          #done
                    "radzimovsky",                      #done
                    "goldsmith",                        #done
                    "dubowsky-freudenstein",            #done
                    "esdu-78035",                       #done
                    "liu et al",                        #done
                    "lankarani-nikravesh cylinder",     #done
                    "pereira et al"]                    #NOK
        return _types

    def _type_check(self):
        """

        :return:
        """
        if self._type.lower() in self._types:
            if self._type.lower() == "pereira et al":
                self._set_parameters_Pereira_etal()
        else:
            raise ValueError, "Cylindrical contact model type not correct! Define contact model type."

    def _evaluate_R(self, R_i=None, R_j=None, dR=None):
        """

        :return:
        """
        if R_i is not None:
            self.R_i = R_i

        if R_j is not None:
            self.R_j = R_j

        if dR is not None:
            self.dR = dR

        if self.R_i is not None and self.R_j is not None and self.dR is not None:
            self.R = (self.R_i * self.R_j) / self.dR
        else:
            self.R = None

    def _evaluate_dR(self, R_i=None, R_j=None):
        """

        :return:
        """
        if R_i is not None:
            self.R_i = R_i
        elif self._parent is not None:
            if self._parent._contact_type.lower() == "revolute clearance joint":
                self.R_i = self._parent.R0_i
            if self._parent._contact_type.lower() == "pin-slot clearance joint linear":
                self.R_i = self._parent.h0_iP / 2.

        if R_j is not None:
            self.R_j = R_j
        elif self._parent is not None:
            if self._parent._contact_type.lower() in ["revolute clearance joint", "pin-slot clearance joint linear"]:
                self.R_j = self._parent.R0_j

        if hasattr(self, "R_i") and hasattr(self, "R_j"):
            if self.R_i is not None and self.R_j is not None:
                if self.contact_location == "internal":
                    self.dR = self.R_i - self.R_j
                if self.contact_location == "external":
                    self.dR = self.R_i + self.R_j
        else:
            self.dR = None

    def _evaluate_E(self, E_list=None, ni_list=None):
        """
        Function evaluates a equivalent modulus of elasticity of the contact
        :return:
        """
        self.h_list = []

        if E_list is not None and ni_list is not None:
            #   overwrite object attributes by function input from user
            self.E_list = E_list
            self.ni_list = ni_list

            for _E, _ni in zip(self.E_list, self.ni_list):
                _h = self._h(_E, _ni)

                self.h_list.append(_h)

            self.h_i, self.h_j = self.h_list

        elif hasattr(self._parent, "body_id_list"):
            self.h_i, self.h_j = self._evaluate_h()

        if self.h_list == []:
            self.E = None
        else:
            self.E = (sum(self.h_list))**(-1)

    def _set_parameters_Pereira_etal(self, dR=None):
        """

        :return:
        """
        if dR is not None:
            self.dR = dR

        if self.dR is not None:
            if self.contact_location == "internal":
                self.a = 0.965
                self.b = 0.0965
                Y = self._evaluate_Y()
                self.n = (Y * self.dR)**(-0.005)

            if self.contact_location == "external":
                self.a = 0.39
                self.b = 0.85
                self.n = 1.094

    def _evaluate_Y(self):
        """
        Limits are in mm.
        :return:
        """
        if 0.005 <= self.dR*1E+3 <= 0.34954:    #   dR from m to mm
            # print "self.dR =", self.dR, self.dR * 1E+3
            # print "np.log(self.dR) =", np.log(self.dR), np.log(self.dR * 1E+6)
            Y = 1.51*(np.log((self.dR * 1E+3 * 1E+3)))**(-0.151)
            # print "Y =", Y
        elif 0.34954 <= self.dR*1E+3 <= 10:
            Y = 0.0151*self.dR + 1.151
        else:
            raise ValueError, "Value for attribute dR not in expected range!"

        return Y

    def contact_force(self, delta, dq_n, dq0_n=None):#dq_t, n, t
        """
        Function calculates normal contact force and tangent friction force based on normal contact force
        :param delta:
        :param dq_n:
        :param dq_t:
        :param n:
        :param t:
        :return:
        """
        if self._type in ContactModel._supported_types():
            Fn = self._contact_force(delta, dq_n)
        elif self._type in self._types:
            Fn = self._contact_force_cylinder(delta, dq_n)

            if self._hysteresis_damping_model is not None:
                if self._hysteresis_damping_model in self._hysteresis_damping_models:
                    damping_factor = self._evaluate_hysteresis_damping()
                    #   damping factor is added to contact force
                    Fn = Fn * (1 + damping_factor)

        else:
            raise ValueError, "Contact model type not defined! Define contact model type. Contact model is %s"%self._type
        
        return Fn

    def _evaluate_hysteresis_damping(self):
        """

        :return:
        """
        if self._hysteresis_damping_model == "hunt-crossley":
            damping_factor = self._hysteresis_damping_factor_HC()
        elif self._hysteresis_damping_model == "herbert-mcwhannell":
            damping_factor = self._hysteresis_damping_factor_HM()
        elif self._hysteresis_damping_model == "lee-wang":
            damping_factor = self._hysteresis_damping_factor_LW()
        elif self._hysteresis_damping_model == "lankarani-nikravesh":
            damping_factor = self._hysteresis_damping_factor_LN()
        elif self._hysteresis_damping_model == "gonthier et al":
            damping_factor = self._hysteresis_damping_factor_Gonthier_etal()
        elif self._hysteresis_damping_model == "zhiying-qishao":
            damping_factor = self._hysteresis_damping_factor_ZQ()
        elif self._hysteresis_damping_model == "flores et al":
            damping_factor = self._hysteresis_damping_factor_Flores_etal()
        elif self._hysteresis_damping_model == "hu-guo":
            damping_factor =  self._hysteresis_damping_factor_HG()
        else:
            print "Not supported!"

        return damping_factor

    def _contact_force_cylinder(self, delta, dq_n):

        delta = abs(delta)
        if self._type == "hertz cylinder":
            Fn = self._Fn_Hertz(delta)
        elif self._type == "johnson":
            Fn = self._Fn_Johnson(delta)
        elif self._type == "radzimovsky":
            Fn = self._Fn_Radzimovsky(delta)
        elif self._type.lower() == "goldsmith":
            Fn = self._Fn_Goldsmith(delta)
        elif self._type.lower() == "dubowsky-freudenstein":
            Fn = self._Fn_DubowskyFreudenstein(delta)
        elif self._type.lower() == "esdu-78035":
            Fn = self._Fn_ESDU_78035(delta)
        elif self._type.lower() == "persson":
            Fn = self._Fn_Persson(delta)
        elif self._type.lower() == "liu et al":
            Fn = self._Fn_Liu_etal(delta)
        elif self._type.lower() == "lankarani-nikravesh cylinder":
            Fn = self._Fn_LankaraniNikraveshCylinder(delta)
        elif self._type.lower() == "pereira et al":
            Fn = self._Fn_Pereira_etal(delta, dq_n)
        else:
            raise ValueError, "Type not corrent. Type is %s"%self._type

        return Fn

    def __evaluate_EdRFn(self, Fn):
        """
        Evaluate expression E*dR / Fn
        :return: EdRFn
        """
        EdRFn = (self.E * self.dR)/abs(Fn)
        return EdRFn

    def __evaluaate_hihjL(self):
        """
        Evaluate expression (self.h_i + self.h_j) / self.L
        :return:
        """
        return (self.h_i + self.h_j) / self.L

    def _Fn_Hertz(self, delta):
        """
        Hertz contact force model for contact between two cylinders
        :param delta:
        :return:
        """
        # print "__Fn_Hertz()"
        #   initial approximation value of contact force
        self._Fn0 = self._Fn_last

        #   contact force per unit length
        _Fn = self._Fn_last = scipy.optimize.newton(self.__evaluate_Fn_Hertz, self._Fn0, args=(delta, ), tol=self._tol, maxiter=self._maxiter)
        return _Fn * self.L

    def _evaluate_delta_Hertz(self, Fn):
        """

        :param Fn:  normal contact force
        :return:
        """
        # print "Fn (_evaluate_delta_Hertz) =", Fn
        EdRFn = self.__evaluate_EdRFn(Fn)
        # print "EdRFn =", EdRFn
        delta = ((2. * abs(Fn))/(self.E * np.pi)) * (np.log(np.pi * EdRFn) - 1.)
        return delta

    def __evaluate_Fn_Hertz(self, Fn, delta):
        """
        Hertz contact force model as implicit function
        :return:
        """
        # print "__evaluate_Fn_Hertz()"
        # self.EdRFn = self.__evaluate_EdRFn(Fn)

        f = self._evaluate_delta_Hertz(Fn) - abs(delta)
        # print "f =", f
        return f

    def __evaluate_dFn_Hertz(self, Fn, delta):
        """

        :return:
        """
        EdRFn = self.__evaluate_EdRFn(Fn)
        # print "EdRFn =", EdRFn
        df = (2. / (self.E * np.pi)) * (np.log(EdRFn * np.pi) - 2.)
        # print "df =", df
        return df

    def _Fn_Johnson(self, delta):
        """

        :param delta:
        :return:
        """
        #   initial approximation value of contact force
        self._Fn0 = self._Fn_last

        #   contact force per unit length
        _Fn = self._Fn_last = scipy.optimize.newton(self.__evaluate_Fn_Johnson, self._Fn0, args=(delta, ), tol=self._tol, maxiter=self._maxiter)
        return _Fn * self.L

    def _evaluate_delta_Johnson(self, Fn):
        """

        :param Fn:
        :return:
        """
        EdRFn = self.__evaluate_EdRFn(Fn)

        delta = (Fn / (self.E * np.pi)) * (np.log(4 * np.pi * EdRFn) - 1) #( * (np.log(4. * self.EdRFn) - 1))
        return delta

    def __evaluate_Fn_Johnson(self, Fn, delta):
        """

        :param Fn:
        :param delta:
        :return:
        """
        f = self._evaluate_delta_Johnson(Fn) - delta

        return f

    def __evaluate_dFn_Johnson(self, Fn, delta):
        """

        :param Fn:
        :param delta:
        :return:
        """
        EdRFn = self.__evaluate_EdRFn(Fn)

        df = (-2 + np.log(4 * EdRFn * np.pi)) / (self.E * np.pi)
        return df

    def _Fn_Radzimovsky(self, delta):
        """

        :param delta:
        :return:
        """
        #   initial approximation value of contact force
        self._Fn0 = self._Fn_last
        #   contact force per unit length
        _Fn = self._Fn_last = scipy.optimize.newton(self.__evaluate_Fn_Radzimovsky, self._Fn0, fprime=self.__evaluate_dFn_Radzimovsky, args=(delta, ), tol=self._tol, maxiter=self._maxiter)
        return _Fn * self.L

    def _evaluate_b(self, Fn):
        """

        :return:
        """
        b = (8/5.) * (abs(Fn) * self.R / self.E)**(1/2.)
        return b

    def _evaluate_delta_Radzimovsky(self, Fn):
        """

        :param delta:
        :return:
        """
        b = self._evaluate_b(Fn)

        delta = (abs(Fn) / (self.E * np.pi)) * ((2/3.) + np.log(4*self.R_i/b) + np.log(4*self.R_j/b))
        return delta

    def __evaluate_Fn_Radzimovsky(self, Fn, delta):
        """

        :param Fn:
        :param delta:
        :return:
        """
        f = self._evaluate_delta_Radzimovsky(Fn) - abs(delta)
        return f

    def __evaluate_dFn_Radzimovsky(self, Fn, delta):
        """

        :param Fn:
        :param delta:
        :return:
        """
        b = self._evaluate_b(Fn)

        df = (1/(self.E * np.pi)) * ((2/3.) + np.log(4 * self.R_i / b) + np.log(4 * self.R_j / b))
        return df

    def _Fn_Goldsmith(self, delta):
        """

        :param delta:
        :return:
        """
        #   initial approximation value of contact force
        self._Fn0 = self._Fn_last

        _Fn = self._Fn_last = scipy.optimize.newton(self.__evaluate_Fn_Goldsmith, self._Fn0, fprime=self.__evaluate_dFn_Goldsmith, args=(delta, ), tol=self._tol, maxiter=self._maxiter)
        return _Fn

    def _evaluate_delta_Goldsmith(self, Fn):
        """

        :param Fn:
        :return:
        """
        self.m = 1
        delta = abs(Fn) * ((self.h_i + self.h_j)/self.L) * (np.log(self.L**self.m / (abs(Fn) * self.R * (self.h_i + self.h_j))) + 1)
        return delta

    def __evaluate_Fn_Goldsmith(self, Fn, delta):
        """s

        :param delta:
        :return:
        """
        f = self._evaluate_delta_Goldsmith(Fn) - abs(delta)
        return f

    def __evaluate_dFn_Goldsmith(self, Fn, delta):
        """

        :param Fn:
        :param delta:
        :return:
        """
        __LmPiFnRhihj = (self.L**self.m)/(abs(Fn) * self.R * (self.h_i + self.h_j))

        df = ((self.h_i + self.h_j) / self.L) * np.log(__LmPiFnRhihj)
        return df

    def _Fn_DubowskyFreudenstein(self, delta):
        """
        Dubowsky-Freudenstein contact force model for pin inside a cylinder.
        Function of contact force in implicit function of Fn, delta and iterative technique is used to solve for Fn
        Body i - pin
        Body j - cylinder
        :param delta:
        :return:
        """
        #   initial approximation value of contact force
        self._Fn0 = self._Fn_last
        # print "--------------------------------"
        # print "delta =", abs(delta)
        # print "F0 =", self._Fn0
        # print "initial approximation =", self._Fn0
        _Fn = self._Fn_last = scipy.optimize.newton(self.__evaluate_Fn_DubowskyFreudenstein, self._Fn0, fprime=self.__evaluate_dFn_DubowskyFreudenstein, args=(delta, ), tol=self._tol, maxiter=self._maxiter)
        # print "F =", _Fn
        return _Fn

    def _evaluate_delta_DubowskyFreudenstein(self, Fn):
        """

        :param Fn:
        :return:
        """
        self._f_log_DF = self._evaluate_f_log_DF(Fn)

        delta = (abs(Fn) * ((self.h_i + self.h_j)/self.L) * (np.log(self._f_log_DF) + 1))
        return delta

    def _evaluate_f_log_DF(self, Fn):
        """

        :return:
        """
        #   part of function written separately for clearer view of code
        # print "Fn =", Fn, "self.R_i =", self.R_i, "self.R_j =", self.R_j, "self.h_i =", self.h_i, "self.h_j =", self.h_j
        _f_log_DF = (self.L**self.m * (self.R_i - self.R_j))/(abs(Fn) * self.R_i * self.R_j * (self.h_i + self.h_j))
        return _f_log_DF

    def __evaluate_Fn_DubowskyFreudenstein(self, Fn, delta):
        """
        Dubowsky-Freudenstein contact force model as implicit function
        :return:
        """
        # print "F0 =", Fn
        f = self._evaluate_delta_DubowskyFreudenstein(Fn) - abs(delta)
        # print "f =", f
        return f

    def __evaluate_dFn_DubowskyFreudenstein(self, Fn, delta):
        """
        Analytical derivative of function __evaluate_Fn_DF(), derivative of main function is calculated first in
        method scipy.optimize.newton()
        :param Fn:
        :param delta:
        :return:
        """
        #   part of function written separately for clearer view of code
        self._f_log_DF = self._evaluate_f_log_DF(Fn)

        #   analytical derivative of function
        df = ((self.h_i + self.h_j)/self.L) * (np.log(self._f_log_DF) - 1)
        return df

    def _Fn_ESDU_78035(self, delta):
        """

        Body i - pin
        Body j - cylinder
        :param delta:
        :return:
        """
        #   initial approximation value of contact force
        self._Fn0 = self._Fn_last

        _Fn = self._Fn_last = scipy.optimize.newton(self.__evaluate_Fn_ESDU, self._Fn0, fprime=self.__evaluate_dFn_ESDU, args=(delta, ), tol=self._tol, maxiter=self._maxiter)
        return _Fn

    def _evaluate_delta_ESDU78035(self, Fn):
        """

        :param Fn:
        :return:
        """
        hihjL = self.__evaluaate_hihjL()

        delta = abs(Fn) * np.pi * hihjL * (np.log(4 * (self.R_i - self.R_j) / (abs(Fn) * hihjL)) + 1)
        return delta

    def __evaluate_Fn_ESDU(self, Fn, delta):
        """
        ESDU-78035 contact force model as implicit function
        :return:
        """
        f = self._evaluate_delta_ESDU78035(Fn) - abs(delta)
        return f

    def __evaluate_dFn_ESDU(self, Fn, delta):
        """
        Analytical derivative of function __evaluate_Fn_ESDU(), derivative of main function is calculated first in
        method scipy.optimize.newton()
        :param Fn:
        :param delta:
        :return:
        """
        #   part of function written separately for clearer view of code
        # self.__f_log_ESDU = (4 * self.L * ())/(abs(Fn) * self.R_i * self.R_j * (self.h_i + self.h_j))
        hihjL = self.__evaluaate_hihjL()
        #   analytical derivative of function
        df = (np.pi * (self.h_i + self.h_j)/self.L) * np.log((4 * (self.R_i - self.R_j))/(abs(Fn) * hihjL))
        return df

    def _Fn_Persson(self, delta):
        """
        Not working
        :param delta:
        :return:
        """
        eps = np.arccos(self.dR / (self.dR + abs(delta)))
        b = np.arctan(eps/2.)
        _Fn = (self.E * self.dR * (b**2 + 1) * b**2) / (2 - (np.log(b**2 + 1) + 2*(b**4)))

        return _Fn

    def _Fn_Liu_etal(self, delta):
        """

        :param delta:
        :return:
        """
        #   contact force per unit length
        _Fn = (1/2.) * np.pi * abs(delta) * self.E * (abs(delta) / (2 * (self.dR + abs(delta))))**0.5
        return _Fn * self.L

    def _Fn_LankaraniNikraveshCylinder(self, delta):
        """

        :param delta:
        :return:
        """
        _Fn = (4./3.) * (self.E * self.R**0.5) * (abs(delta)**self.n)
        return _Fn

    def _Fn_Pereira_etal(self, delta, dq_n):
        """

        :param delta:
        :param dq_n:
        :return:
        """
        _Fn = (((self.a * (self.dR) + self.b) * (self.L) * (self.E)) / (self.dR)) * (abs(delta)**self.n) * (1 + ((3. / 4.) * (1 - self.c_r**2) * (dq_n / self._dq0_n))) #
        # print "delta =", delta, "Fn =", _Fn
        return _Fn * self.L

if __name__ == "__main__":

    _type = "kelvin-voigt" #hertz cylinder, johnson, radzimovsky, goldsmith, dubowsky-freudenstein, esdu-78035, persson, liu et al, lankarani-nikravesh, pereira et al
    # print _type
    model = ContactModelCylinder(_type)
    #   hole/slot
    model.R_i = 2.5E-3
    #   pin
    model.R_j = 2.45E-3
    #   contact length
    model.L = 1
    #   exponent
    model.m = 3
    #   max iter
    model._maxiter = 10

    #   material properties: aluminium
    model._evaluate_E([69E+9, 69E+9], [0.35, 0.35])

    #   additional geometrical properties
    model._evaluate_dR()
    model._evaluate_R()
    # model._set_parameters_Pereira_etal()
    pprint(vars(model))

    model.c_r = .5
    _delta = -1E-5
    _dqn0 = 0.
    model._dq0_n = -1E-2
    # print model.contact_force(_delta, _dqn0)

    mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    mpl.rcParams['text.usetex'] = False #True, False

    fig = plt.figure(num=1, figsize=(6, 4), dpi=100, facecolor='w', edgecolor='k')
    ax = plt.subplot(111, aspect="auto")
    ax.ticklabel_format(style='sci',scilimits=(-4,4), axis='both')
    ax.set_title(model._type)
    plt.xlabel(r'$\delta$ [m]', fontsize=10)
    plt.ylabel(r'$dq_n$ [m/s]', fontsize=10)

    #   vector of input  variable
    delta_iter = np.arange(1E-10, 1E-6, 5E-9)
    Fn_iter = np.zeros(len(delta_iter))

    #   solve for delta
    Fn = np.linspace(1E-6, 1E+4, len(delta_iter))
    delta = np.zeros(len(Fn))
    for i in range(0, len(delta)):
        if model._type == "hertz":
            delta[i] = model._evaluate_delta_Hertz(Fn[i])
        if model._type == "johnson":
            delta[i] = model._evaluate_delta_Johnson(Fn[i])
        if model._type == "radzimovsky":
            delta[i] = model._evaluate_delta_Radzimovsky(Fn[i])
        if model._type == "goldsmith":
            delta[i] = model._evaluate_delta_Goldsmith(Fn[i])
        if model._type == "dubowsky-freudenstein":
            delta[i] = model._evaluate_delta_DubowskyFreudenstein(Fn[i])
        if model._type == "esdu-78035":
            delta[i] = model._evaluate_delta_ESDU78035(Fn[i])

    plt.plot(delta, Fn, "--", label="solve for delta")
    # print "delta =", delta

    #   solve for Fn iteratively
    for i in range(0, len(delta_iter)):
        Fn_iter[i] = abs(model.contact_force(delta_iter[i], _dqn0))
        # print "delta =", delta_iter[i], "Fn_iter =", Fn_iter[i]

    # print "   delta(in),       Fn(out),         Fn(in),          delta(out)"
    # print np.hstack((delta_iter,
    #                 Fn_iter,
    #                 Fn,
    #                 delta))

    plt.plot(delta_iter, Fn_iter, "-", label="solve for Fn-iter")

    plt.legend(loc="best")
    # plt.xlim([-max(delta), max(delta)])
    # plt.ylim([-max(Fn), max(Fn)])
    print "Fn =", model.contact_force(1E-7, _dqn0)
    plt.show()
    #   save data to excel
    filename = "model_"+_type+".xlsx"
    workbook = xlsxwriter.Workbook(filename)
    worksheet = workbook.add_worksheet("data")

    #   delta
    worksheet.write(0, 0, "delta")
    worksheet.write_column(1, 0, delta_iter)

    #   Fn
    worksheet.write(0, 1, "Fn")
    worksheet.write_column(1, 1, Fn_iter)

    workbook.close()



