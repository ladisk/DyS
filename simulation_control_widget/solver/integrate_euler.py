"""
Created on 10. maj 2016

@author: luka.skrinjar
"""
import copy
import time
import datetime
import numpy as np


from simulation_control_widget.solver.dynamic_analysis import DynamicAnalysis
from DAE_fun import DAEfun


class IntegrateEuler(DynamicAnalysis):
    """
    classdocs
    """

    def __init__(self, MBD_system, parent=None):
        """
        Constructor
        """
        super(IntegrateEuler, self).__init__(MBD_system, parent)

    def evaluate_error(self, t, h):
        """

        :return:
        """
        e = np.exp(t) - (1. + h)**(t/h)

        return e

    def solve(self):
        """
        Euler method for time integration
        :param t_0:
        :param t_n:
        :param q0:
        :param h:
        """
        #   solution containers
        self._solution_containers()

        #    copy MBD object
        self._MBD_system = copy.copy(self.MBD_system)
        #    ode fun object
        self.DAE_fun = DAEfun(self._MBD_system, parent=self)

        self.start_solve()

        self.start_time_simulation_info_in_sec_UTC = time.time()
        self.start_time_simulation_info = datetime.datetime.fromtimestamp(self.start_time_simulation_info_in_sec_UTC).strftime("%H:%M:%S %d.%m.%Y")
        print "Simulation started: ", self.start_time_simulation_info

        #    create array of initial conditions for differential equations
        #   get initial conditions
        q0 = self._MBD_system.evaluate_q0()
        self.FLAG = 1

        self.simulation_id = 0
        #    0 - no contact
        #    +1 - contact detected
        #    -1 - contact already happened - reduce integration step
        self.step = 0
        self.h_contact = self._MBD_system.Hmin
        self._append_to_file_step = 0

        #    constant time step vector
        self.t_constant_time_step = np.arange(0, self._MBD_system.t_n + self._MBD_system.Hmax, self._MBD_system.Hmax)

        self.t_level = 0
        self.t_0 = self._dt = self._MBD_system.t_0
        self.t_n = self._MBD_system.t_n
        t = self.t_0
        w = q0

        self.Hmax = h = self._MBD_system.Hmax
        self.R = 0

        while self.FLAG == 1:
            if self.stopped:
                self.refresh(t=t, q=w)
                self.stop_solve()
                self.FLAG = 0

            if t >= self.t_n:
                self.FLAG = 0
                self.finished = True
                self.refresh(t=t, q=w)

            self.t = t = t + h
            self.h = h
            # print "-------------------------------------"
            # print self.t, t, h

            #   error control
            self.absError = self.evaluate_error(self.t, self.h)
            if self.absError <= self.absTol:
                w = w + h * self.DAE_fun.evaluate_dq(t, w)

                #    solve contacts
                self.FLAG_contact = self._evaluate_contacts(t, w)

                #   evaluate mechanical energy of system
                self._energy_t = self._mechanical_energy(w)
                self._energy_delta = self._energy_t - self._solution_data._energy_solution_container[-1]
            else:
                self.t = t = t - h
                self.h = h / 2.

            #   track data
            h, t, w = self._track_data(h, t, w)

            #   simulation informations
            self._info(t, w)

            #   check time step
            h, self._error = self.check_time_step(t, h, w)

        if self.finished or self.failed or self.stopped:
            self.finished_solve()

        #    save data to file
        # if self.MBD_system._solution_save_options.lower() != "discard":
        self.write_simulation_solution_to_file()

    def check_time_step(self, t, h, w):
        """

        :param t:
        :param h:
        :param w:
        :return:
        """
        error = False

        if self.absError <= self.absTol:
            h = 2. * h

        else:
            error = True

        return h, error

