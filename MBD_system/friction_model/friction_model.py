"""
Created on 9. jul. 2014

@author: lskrinjar
"""
from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib import rcParams


from MBD_system.fix_string import fix_string


class FrictionModel(object):
    """
    classdocs
    """
    def __init__(self, _type, coef_of_friction_kinematic=None, coef_of_friction_dynamic=None, coef_of_friction_static=None, v0=0, v1=0, v_sigma=0, delta_v=1, gamma=1, properties_dict={}, parent=None):
        """
        Constructor of friction model class
        supported input for parameter
        :param _type:   supported types of friction models
                        1 Coulomb
                        2 Viscous - check
                        3 Stribeck
                        4 Karnopp
                        5 Threlfall
                        6 Dahl
                        7 Ambrosio
                        8 Sjo
                        8 LuGre - todo
                        9 Leuven - todo
                        10 Bristle - todo
                        11 Reset integrator - todo
                        12 Karnopp - todo
                        13 Bliman-Sorine - todo
        :param coef_of_friction_kinematic:
        :param coef_of_friction_static:
        :param coef_of_friction_dynamic:
        :param v0:
        """
        #   parent
        self._parent = parent

        #   type
        self._type = _type
        self._types = ["coulomb",
                       "viscous",
                       "karnopp",
                       "stribeck",
                       "threlfall",
                       "dahl",
                       "ambrosio",
                       "lugre",
                       "leuven",
                       "bristle",
                       "reset integrator",
                       "karnopp",
                       "bliman-sorine"]

        #   sub string
        self._substring = "friction_model"

        #   friction coefficients
        self.coef_of_friction_static = coef_of_friction_static
        self.coef_of_friction_kinematic = coef_of_friction_kinematic
        self.coef_of_friction_dynamic = coef_of_friction_dynamic

        #   stiction friction model

        #   viscous friction model
        #   coefficient of viscous friction
        #   stribeck velocity - sliding speed coefficient
        self.v_sigma = v_sigma
        #   gradient of friction decay
        self.gamma = gamma
        #   non-linear dependence on velocity
        self.delta_v = delta_v

        #   slip velocity values
        self.v0 = v0
        self.v1 = v1
        if self._type.lower() == "ambrosio":
            if v1 < v0:
                self.v0 = v1
                self.v1 = v0

        #   tolerance
        self._dq_t_TOL = 1E-3

        self._add_aditional_parameters(properties_dict)

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

    def setType(self, type):
        """

        :param type:
        :return:
        """
        self._type = fix_string(type.title())

    def friction_force(self, Fn, dq_t):
        """
        Function evaluates tangential friction force

        :param Fn:      normal force value (type:float)
        :param dq_t:    relative tangential velocity (type:float)
        :return Ft:     tangential force value (type:float)
        """
        #   get absolute value of contact force
        Fn = abs(Fn)

        if self._type.lower() == "ideal":
            Ft = 0
        elif self._type.lower() == "coulomb":
            Ft = self.__Ft_coulomb(Fn, dq_t)
        elif self._type.lower() == "viscous":
            Ft = self.__Ft_viscous(Fn, dq_t)
        elif self._type.lower() == "stribeck":
            Ft = self.__Ft_stribeck(Fn, dq_t)
        elif self._type.lower() == "karnopp":
            Ft = self.__Ft_karnopp(Fn, dq_t)
        elif self._type.lower() == "threlfall":
            Ft = self.__Ft_threlfall(Fn, dq_t)
        elif self._type.lower() == "dahl":
            Ft = self.__Ft_dahl(Fn, dq_t)
        elif self._type.lower() == "ambrosio":
            Ft = self.__Ft_ambrosio(Fn, dq_t)
        elif self._type.lower() == "sjo":
            Ft = self.__Ft_sjo(Fn, dq_t)
        else:
            raise ValueError, "Contact model type not defined! Define contact model type. User defined type is:%s" %(self._type)

        return Ft

    def __Ft_coulomb(self, Fn, dq_t):
        """
        Coulomb friction model
        :param Fn:      normal contact force (type float, non-negative)
        :param _dq_t:
        :return:
        """
        if abs(dq_t) <= self._dq_t_TOL:
            Ft = 0.
        else:
            Ft = self.coef_of_friction_kinematic * Fn * np.sign(dq_t)
        return Ft

    def __Ft_viscous(self, Fn, dq_t):
        """
        Viscous friction model with coulomb friction
        :param Fn:
        :param _dq_t:
        :return:
        """
        Fc = self.__Ft_coulomb(Fn, dq_t)

        Ft = Fc + np.sign(dq_t) * self.v_sigma * abs(dq_t)
        return Ft

    def __Ft_stribeck(self, Fn, dq_t):
        """
        Stribeck friction model
        :param Fn:
        :param _dq_t:
        :return:
        """
        Fc = self.__Ft_coulomb(Fn, dq_t)
        Fs = self.coef_of_friction_static * Fn * np.sign(dq_t)
        Ft = Fc + (Fs - Fc) * np.exp(-abs(dq_t/self.v_sigma)**self.delta_v) + self.v_sigma * dq_t
        return Ft

    def __Ft_threlfall(self, Fn, dq_t):
        """

        :param Fn:
        :param dg_t:
        :return:
        """
        if abs(dq_t) <= self.v0:
            coef_of_friction = self.coef_of_friction_kinematic * (1 - np.exp(-3 * abs(dq_t) / self.v0))
        else:
            coef_of_friction = 0.95 * self.coef_of_friction_kinematic

        Ft = coef_of_friction * Fn * np.sign(dq_t)
        return Ft

    def __Ft_dahl(self, Fn, dq_t):
        """

        :param Fn:
        :param dq_t:
        :return:
        """
        if abs(dq_t) <= self.v0:
            coef_of_friction = self.coef_of_friction_kinematic * abs(dq_t) / self.v0
        else:
            coef_of_friction = self.coef_of_friction_kinematic

        Ft = coef_of_friction * Fn * np.sign(dq_t)
        return Ft

    def __Ft_ambrosio(self, Fn, dq_t):
        """
        Ambrosio friction model
        :param Fn:
        :param dq_t:
        :return:
        """
        if abs(dq_t) <= self.v0:
            coef_of_friction = 0
        elif self.v0 < abs(dq_t) <= self.v1:
            coef_of_friction = self.coef_of_friction_kinematic * ((abs(dq_t) - self.v0) / (self.v1 - self.v0))
        else:
            coef_of_friction = self.coef_of_friction_kinematic

        Ft = coef_of_friction * Fn * np.sign(dq_t)
        return Ft

    def __Ft_sjo(self, Fn, dq_t):
        """

        :param Fn:
        :param dq_t:
        :return:
        """
        coef_of_friction = self.coef_of_friction_kinematic * np.arctan(self.gamma * abs(dq_t) * np.pi / (2. * self.coef_of_friction_kinematic)) * (2. / np.pi)

        Ft = coef_of_friction * Fn * np.sign(dq_t)
        return Ft

    def __Ft_karnopp(self, Fn, dq_t):
        """

        :param Fn:
        :param dq_t:
        :return:
        """
        if abs(dq_t) > self.v1:
            Ft = self.__Ft_coulomb(Fn, dq_t)
        else:
            #   here has to be implemented to access to Qe of a body and tangemt
            Ft = 0.

        return Ft


if __name__ == "__main__":
    types = ["Coulomb",
            "Viscous",
            "Karnopp",
            "Stribeck",
            "Threlfall",
            "Dahl",
            "Ambrosio",
            "Sjo"]

    coef_of_friction_kinematic = 0.5
    coef_of_friction_dynamic = 0.2
    coef_of_friction_static = 0.6

    v0=10
    v1=4
    v_sigma=.002
    delta_v=2

    models = []

    Fn = 1
    dv = 0.001
    v = np.arange(-40, 40+dv, dv)
    gamma = .1

    for _type in types:
        model = FrictionModel(_type,
                              coef_of_friction_kinematic=coef_of_friction_kinematic,
                              coef_of_friction_dynamic=coef_of_friction_dynamic,
                              coef_of_friction_static=coef_of_friction_static,
                              v0=v0,
                              v1=v1,
                              v_sigma=v_sigma,
                              delta_v=delta_v,
                              gamma=gamma)

        models.append(model)

    sb.set_style("whitegrid")
    fig = plt.figure(num=1, figsize=(6, 4), dpi=100, facecolor='w', edgecolor='k')
    ax = plt.subplot(111, aspect="auto")
    ax.ticklabel_format(style='sci',scilimits=(-4,4), axis='both')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    #   tex settings
    rcParams['text.usetex'] = True
    rc('font', family='serif')

    for i, model in enumerate(models):
        model.Ft_solution_data = []

        for v_i in v:
            Ft = model.friction_force(Fn, v_i)

            model.Ft_solution_data.append(Ft)

        if i in [0, 1, 2, 3, 4, 5, 6, 7]:
            plt.plot(v, np.array(model.Ft_solution_data), label=model._type)


    plt.xlabel(r'\ $v_t$', fontsize=12)
    plt.ylabel(r'\ $\mu$',  fontsize=12)
    plt.ylim([-1., 1.])
    legend = plt.legend(frameon=1, loc='upper left')
    frame = legend.get_frame()
    frame.set_color('white')


    # plt.savefig("c:\Users\lskrinjar\Documents\FS_3_pd\podiplomski seminar 2\PorociloPodiplomskiSeminar2\teoreticneOsnove\modeliSileTrenja\modeliSileTrenja_primerjava.pdf")
    filename = "modeliSileTrenja_primerjava.pdf"
    plt.savefig(filename)

    # plt.show()

