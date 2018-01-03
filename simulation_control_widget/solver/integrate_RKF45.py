"""
Created on 10. maj 2016

@author: luka.skrinjar
"""
import copy
import time
import datetime
import numpy as np
import logging
from pprint import pprint
from PyQt4 import QtCore, QtGui


from simulation_control_widget.solver.dynamic_analysis import DynamicAnalysis
from DAE_fun import DAEfun
from ODE_fun import ODE_fun


class IntegrateRKF45(DynamicAnalysis):
    """
    classdocs
    """

    def __init__(self, MBD_system, parent=None):
        """
        Constructor
        """
        super(IntegrateRKF45, self).__init__(MBD_system, parent)

    def evaluate_absError(self, w, w_new):
        """

        :return:
        """
        # e = (1. / self.h) * np.linalg.norm(((1./360.)*K1 - (128./4275.)*K3 - (2197./75240.)*K4 + (1./50.)*K5 + (2./55.)*K6), ord=2)

        # w_new = self._evaluate_w_new(w, K1, K2, K3, K4, K5, K6)

        e = self.evaluate_error_HHT_I3(w, w_new)

        # if e < 1E-11:
        #     e = 1E-11

        return e

    def evaluate_error_HHT_I3(self, w, w_new):
        """

        :return:
        """
        q_new = self.get_q(w_new)

        self.evaluate_q_max(w_new)

        x = q_new / self.q_max

        #   error
        e = (self.h**2 / np.sqrt(self.q_n)) * np.linalg.norm(x, ord=2)#(self.h**2 / np.sqrt(self.q_n)) * np.linalg.norm(x, ord=2)
        return e

    def evaluate_error_RKF45(self, K1, K2, K3, K4, K5, K6):
        """

        :return:
        """
        d = (1./self.h) * np.abs(np.array((1./360.)*K1 - (128./4275.)*K3 - (2197./75240.)*K4 + (1./50.)*K5 + (2./55.)*K6))

        e = np.max(np.abs(self.get_q(d)))
        return e

    def _evaluate_w_new(self, w, K1, K2, K3, K4, K5, K6):
        """

        :return:
        """
        w_new = w + (25. / 216.) * K1 + (1408. / 2565.) * K3 + (2197. / 4104.) * K4 - (1. / 5.) * K5

        return w_new

    def check_error(self):
        """

        :return:
        """
        if self.absError is np.nan:
            self.simulation_error.setError("Evaluated error is NaN!")
            return True

        if self.absError > self.absTol:
            return True
        else:
            return False

    def check_time_step(self, t, h, w, R=None):
        """

        :param t:
        :param h:
        :param w:
        :return:
        """
        if self.FLAG_contact == 0:
            delta = 0.84 * (self.absTol / np.max(R)) ** (1./4.)

            if delta <= 0.1:
                h = 0.1 * h

            elif delta >= 4.:
                h = 4. * h

            else:
                h = delta * h

            if h > self.Hmax:
                h = self.Hmax

            if t >= self.t_n:
                self.FLAG = 0

            elif t + h > self.t_n:
                h = self.t_n - t

            elif h < self.Hmin:
                print "h =", h
                print "self.Hmin =", self.Hmin
                print "self.absTol =", self.absTol
                print "self.absError =", self.absError
                # h = self.Hmin
                self.simulation_error.setWarning("Hmin exceeded!")
                # error = True
                # h = self.Hmin
                # print "Hmin exceeded!"
                # self.FLAG = 0

        if self.FLAG_contact == 1:
            h = self.Hcontact

        if self.FLAG_contact == -1:
            h = 0.5 * h

        if h < self.Hmin and t < self.t_n - self.Hmin:
            print "h =", h
            print "self.Hmin =", self.Hmin
            print "self.absTol =", self.absTol
            print "self.absError =", self.absError
            print "self.contact_status_list =", self.contact_status_list
            self.simulation_error.setError("Hmin exceeded! FLAG contact is: " + str(self.FLAG_contact), q=w)
            # error = True
            # h = self.Hmin

        return h, self.simulation_error.error

    def solve(self, q0=[], t_n=None, absTol=1E-9):
        """
        Solves system of ode with order 5 method with runge-kutta algorithm
        """
        #   redefine solution containers
        self._solution_containers()

        #    copy MBD object
        if hasattr(self.MBD_system, "setSolver"):
            self.MBD_system.setSolver(self)
            self._MBD_system = self.MBD_system
            #   set solver
            self._MBD_system.setSolver(self)

        #    dae fun object
        if self.DAE_fun is None:
            self.DAE_fun = DAEfun(self._MBD_system, parent=self)

        if q0:
            self.q0 = q0

        self.start_solve()

        #   create array of initial conditions for differential equations
        #   get initial conditions
        if hasattr(self._MBD_system, "q0"):
            self._MBD_system.q0 = self.q0

        if t_n is not None:
            self.t_n = t_n

        #   evaluate stiffness matrix
        if hasattr(self.DAE_fun, "evaluate_K"):
            self.DAE_fun.evaluate_K(self.q0)

        self.simulation_id = 0
        #    0 - no contact
        #    +1 - contact detected
        #    -1 - contact already happened - reduce integration step
        self.step = 0
        if hasattr(self._MBD_system, "Hmin"):
            self.h_contact = self._MBD_system.Hmin

        self._append_to_file_step = 0

        self.t_0 = self._dt = 0
        self.FLAG = 1

        if hasattr(self._MBD_system, "t_n"):
            if self._MBD_system.t_n is not None:
                self.t_n = self._MBD_system.t_n

        if hasattr(self._MBD_system, "stepsNumber"):
            self.stepsNumber = self._MBD_system.stepsNumber

        t = self.t_0
        w = self.q0

        if self.t_n is None and self.stepsNumber is not None:
            self.t_n = np.inf

#         logging.getLogger("DyS_logger").info("Size of MBD system in bytes %s" % sys.getsizeof(self.MBD_system))
        logging.getLogger("DyS_logger").info("Simulation started with initial conditions q0:\n%s" % self.q0)

        h = self.Hmax
        self._t_FLAG1 = self.Hmax

#        print "absTol =", absTol
#        self.update_opengl_widget_every_Nth_step = 1*((t_n - t_0)/Hmax)
#        print "self.update_opengl_widget_every_Nth_step =", self.update_opengl_widget_every_Nth_step

        self.save_updated_screenshot_to_file(self.t_0)
        # np.set_printoptions(precision=20, threshold=1000, edgeitems=True, linewidth=1000, suppress=False, formatter={'float': '{: 10.9e}'.format})

        while self.FLAG == 1 and not self._error and not self.DAE_fun.error:
            # time.sleep(.001)
            # if self.step > 1415:
            #     print "----------------------"
            #     print "step =", self.step, "t =", t, "h =", h#"self.FLAG_contact =", self.FLAG_contact
            # print "step =", self.step, "h =", h, "t =", t
            # print "step =", self.step
            self.h = h
            self.t = t
            # print "self.h =", self.h
            K1 = h * self.DAE_fun.evaluate_dq(t, w)
            K2 = h * self.DAE_fun.evaluate_dq(t + (1. / 4.) * h, w + (1. / 4.) * K1)
            K3 = h * self.DAE_fun.evaluate_dq(t + (3. / 8.) * h, w + (3. / 32.) * K1 + (9. / 32.) * K2)
            K4 = h * self.DAE_fun.evaluate_dq(t + (12. / 13.) * h, w + (1932. / 2197.) * K1 - (7200. / 2197.) * K2 + (7296. / 2197.) * K3)
            K5 = h * self.DAE_fun.evaluate_dq(t + h, w + (439. / 216.) * K1 - 8. * K2 + (3680. / 513.) * K3 - (845. / 4104.) * K4)
            K6 = h * self.DAE_fun.evaluate_dq(t + (1. / 2.) * h, w - (8. / 27.) * K1 + 2. * K2 - (3544. / 2565.) * K3 + (1859. / 4104.) * K4 - (11. / 40.) * K5)

            w_new = self._evaluate_w_new(w, K1, K2, K3, K4, K5, K6)

            # self.absError = (1. / h) * np.linalg.norm((1. / 360.) * K1 - (128. / 4275.) * K3 - (2197. / 75240.) * K4 + (1. / 50.) * K5 + (2. / 55.) * K6)

            # self.absError = self.evaluate_error_HHT_I3(w, w_new, K1, K2, K3, K4, K5, K6)
            self.absError = self.evaluate_absError(w, w_new)

            #self.absError = (1. / h) * np.linalg.norm((1. / 360.) * K1 - (128. / 4275.) * K3 - (2197. / 75240.) * K4 + (1. / 50.) * K5 + (2. / 55.) * K6)
            #self.absError = 0.1*self.absTol
            # self.R = absTol
            # print "E_RKF45 =", self.evaluate_error_RKF45(K1, K2, K3, K4, K5, K6), "E_HHT_I3 =", self.absError
            #    if calculated difference is less the absolute tolerance limit and accept the calculated solution
            # if not self.errorControl:
            #     self.absError = self.absTol

            # if (self.absError <= self.absTol) and (self.FLAG_contact in [0, 1]):
            # print "self.check_error() =", self.check_error()

            #    solve contacts
            self.FLAG_contact = self._evaluate_contacts(t, w_new)

            if not self.check_error():#and (self.FLAG_contact in [0, 1])
                #   next time point
                self.t = t = t + h
                #   value of vector of variables at next time step
                w = w_new

                if self.FLAG_contact in [0, 1]:
                    #   evaluate mechanical energy of system
                    self._mechanical_energy, self._kinetic_energy, self._potential_energy, self._elastic_strain_energy = self._evaluate_mechanical_energy(w)
                    self._energy_delta = self.evaluate_energy_delta(self._mechanical_energy)

                else:
                    self.t, self.step, self._mechanical_energy, self._energy_delta = self._use_last_accapted_solution(t, self.step)

            else:
                self.t, self.step, self._mechanical_energy, self._energy_delta = self._use_last_accapted_solution(t, self.step)

            t = self.t

                # w = w
                # self.step = self.step
                # self._mechanical_energy = 0.
                # self._energy_delta = 0.

            #    track - store simulation state data of different MBD system objects
            # print "h1 =", h
            h, t, w = self._track_data(h, t, w)
            # print "h2 =", h
            #   check time step
            if self.errorControl:
                h, self._error = self.check_time_step(t, h, w, R=self.absError)
            else:
                print UserWarning, "Error control dissabled! Results may be wrong!"
                self._error = False
            # print "h3 =", h
            #    process information of integration process
            if not self.check_error() and (self.FLAG_contact in [0, 1]) and self.update_visualization_widget:
                self._info(t, w)

            #   check if finished
            if (t >= self.t_n) or (self.step == self.stepsNumber):
                self.FLAG = 0
                self.finished = True
                self.refresh(t=t, q=w)

            #    integration finished - end time reached successfully
            if self.finished or self.failed:
                self.refresh(t=t, q=w)
                self.finished_solve()
                self.FLAG = 0

            if self.stopped:
                self.refresh(t=t, q=w)
                self.stop_solve()
                self.FLAG = 0

            if np.isnan(t):
                self._error = True
                self.failed = True
                print "t is NaN! Check Hmax!"

        if self._error or self.failed or self.simulation_error.error:
            print "self._error =", self._error
            print "self.failed =", self.failed
            print "self.DAE_fun.error =", self.DAE_fun.error
            self.simulation_failed()
            print self.simulation_error.info
            print Warning, "Simulation failed!"

        #    save data to file
        if self.write_to_file:
            self.write_simulation_solution_to_file()

    def _use_last_accapted_solution(self, t, step):
        """

        :return:
        """
        energy_t = 0.
        energy_delta = 0.

        return t, step, energy_t, energy_delta


def evaluate_dq(t, y):
    """

    :param t:
    :param y:
    :return:
    """
    dy = y - t**2 + 1.
    return dy


def y_solution(t):
    """

    :param t:
    :return:
    """
    y = (t + 1.)**2 - 0.5 * np.exp(t)
    return y


if __name__ == "__main__":
    t_n = 2.
    y0 = 0.5
    Hmax = 0.25
    Hmin = 0.01
    TOL = 1E-5

    int = IntegrateRKF45(None)
    int.DAE_fun = ODE_fun(None)
    int.DAE_fun.evaluate_dq = evaluate_dq
    int.Hmax = Hmax
    int.Hmin = Hmin
    int.absTol = TOL
    int.update_visualization_widget = False
    int.write_to_file = False
    int.solve(y0, t_n, absTol=TOL)

    print "step \tt_i \t\t\ty_i \t\t\tw_i \t\t\th_i \t\t\tR_i \t\t\t|q_i - w_i|"
    for i, t in enumerate(int._solution_data._t_solution_container):
        y_i = y_solution(t)
        print "%1i \t\t%1.7f \t\t%1.7f \t\t%1.7f \t\t%1.7f \t\t%1.7f \t\t%1.7f" % (i, t, y_i, int._solution_data._q_solution_container[i], int._solution_data._h_solution_container[i], int._solution_data._absError_solution_container[i], abs(y_i - int._solution_data._q_solution_container[i]))

