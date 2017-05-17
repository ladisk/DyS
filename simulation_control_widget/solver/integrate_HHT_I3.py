"""

created by: lskrinjar
date of creation: 30/01/2017
time of creation: 17:47
"""

import copy
import time
import datetime
import numpy as np
import logging
from pprint import pprint
from matplotlib import pyplot as plt
from scipy import optimize


from simulation_control_widget.solver.dynamic_analysis import DynamicAnalysis
from DAE_fun import DAEfun
from ODE_fun import ODE_fun


class IntegrateHHT_I3(DynamicAnalysis):
    """
    classdocs
    """

    def __init__(self, MBD_system, parent=None):
        """
        Constructor
        """
        super(IntegrateHHT_I3, self).__init__(MBD_system, parent)

        #   iterative process properties
        self.FLAG_iter = None
        self.J = None

        #   parameters
        self.gamma = 0.5
        self.beta = 0.25
        self.alpha = 0.

    def evaluate_beta(self, alpha):
        """

        :return:
        """
        beta = ((1. - alpha)**2) / 4.
        return beta

    def evaluate_alpha(self, alpha):
        """
        If alpha is 0, no numerical damping is used
        :param alpha:
        :return:
        """
        if alpha < -0.3:
            alpha = -0.3
        if alpha > 0.:
            alpha = 0.

        return alpha

    def evaluate_gamma(self, alpha):
        """

        :param alpha:
        :return:
        """
        gamma = 0.5 + alpha
        return gamma

    def evaluate_error(self, q, q_new):
        """

        :return:
        """
        self.evaluate_q_max(q_new)

        x = q_new / self.q_max

        #   error
        e = abs(self.beta - (1. / (6. * (1. + self.alpha)))) * (self.h**2 / np.sqrt(self.q_n)) * np.linalg.norm(x, ord=2)

        return e

    def solve(self, q0=[], t_n=None, absTol=1E-9, maxIter=10):
        """

        :return:
        """
        #   redefine solution containers
        self._solution_containers()

        #    copy MBD object
        self.MBD_system.setSolver(self)
        # self._MBD_system = copy.copy(self.MBD_system)
        self._MBD_system = self.MBD_system
        #   set solver
        self._MBD_system.setSolver(self)
        #    ode fun object
        self.DAE_fun = DAEfun(self._MBD_system, parent=self)

        self.start_solve()

        self.beta = self.evaluate_beta(self.alpha)
        self.gamma = self.evaluate_gamma(self.alpha)

        self.t_n = t_n
        if self.t_n is not None and self.stepsNumber is not None:
            self.Hmax = self.t_n / self.stepsNumber

        if q0:
            self.q0 = q0

        if t_n is not None:
            self.t_n = t_n

        self.h = h = self.Hmax
        t = self.t_0
        w = self.q0

        q = self.get_q(w)

        while self.FLAG == 1 and not self._error and not self.DAE_fun.error:
            self.h = h
            K1 = w + 0.5 * h * self.DAE_fun.evaluate_dq(t, w)
            w0 = K1

            i_iter = 0
            self.FLAG_iter = 0
            while self.FLAG_iter == 0:
                w = w0 - ((w0 - (0.5 * h * self.DAE_fun.evaluate_dq(t + h, w0)) - K1) / (1. - 0.5 * h * self.DAE_fun.evaluate_Jacobian(t + h, w0)))

                self.absError = self.evaluate_error(w0, w)

                #optimize.newton(f, x0, fprime=None, args=(y,), tol=1.48e-08, maxiter=50)

                if self.absError < self.absTol:
                    self.FLAG_iter = 1
                else:
                    i_iter += 1
                    w0 = w
                    if i_iter > maxIter:
                        raise ValueError, "The maximum number of iterations exceeded"
                        break

            if (t >= self.t_n) or (self.step == self.stepsNumber):
                self.FLAG = 0
                self.finished = True
                self.refresh(t=t, q=w)
            else:
                t = t + h

            #    track - store simulation state data of different MBD system objects
            h, t, w = self._track_data(h, t, w)

    def evaluate_Jacobian(self, t, q):
        """

        :return:
        """
        if self.J is None:
            self.J = self.DAE_fun.evaluate_Jacobian(t, q)

        return self.J


if __name__ == "__main__":
    int = IntegrateHHT_I3(None)
    y0 = -1.
    h = 0.2
    t_n = 1.
    int.DAE_fun = ODE_fun(None)
    int._active_scene = False
    int.Hmax = 0.2
    int.absTol = 1E-6
    int.solve(y0, t_n, absTol=1E-6, maxIter=10)

    print "int =", int._solution_data
    pprint(vars(int._solution_data))
    t = int._solution_data._t_solution_container

    plt.plot(t, int._solution_data._q_solution_container, label="HHT-I3")

    y = np.zeros_like(t)
    for i in range(0, len(t)):
        y[i] = t[i] - np.exp(-5.*t[i])

    print "Solution"
    for i, (t_i, y_i, y_num_i) in enumerate(zip(t, y, int._solution_data._q_solution_container)):
        print i, t_i, y_i, y_num_i

    plt.plot(t, y, label="analytical sol.")

    plt.legend()
    plt.show()


