"""

created by: lskrinjar
date of creation: 16/06/2016
time of creation: 17:57
"""
import copy
import time
import multiprocessing
import threading

from kinematic_analysis import KinematicAnalysis
from integrate_euler import IntegrateEuler
from integrate_RKF45 import IntegrateRKF45


class SolverThread(threading.Thread):
    """

    """

    def __init__(self, MBD_system, queue, lock):
        threading.Thread.__init__(self)

        #   name
        self._name = multiprocessing.current_process().name

        #   MBD system
        self.MBD_system = MBD_system

        #   queue
        self.queue = queue

        #   lock
        self.lock = lock

        #   solver (analysis)
        self.analysis = self._analysis_type()

        #   open GL properties
        self._update_display_type = "dt"
        self._update_display_types = ["step", "dt"]

        self._dt = self.MBD_system._dt

    def _analysis_type(self):
        """

        :return:
        """
        #   kinematic analysis
        if self.MBD_system.analysis_type == "kinematic":
            analysis = KinematicAnalysis(self.MBD_system, parent=self)

        #   dynamic analysis
        elif self.MBD_system.analysis_type == "dynamic":

            if self.MBD_system.integrationMethod == "euler":
                analysis = IntegrateEuler(self.MBD_system, parent=self)

            elif self.MBD_system.integrationMethod == "RKF45":
                analysis = IntegrateRKF45(self.MBD_system, parent=self)

        else:
            analysis = None
            raise ValueError, "Analysis type not defined!"

        return analysis

    def run(self):
        """

        :return:
        """
        # with self.lock:
        #   update all variables
        self._set_variables()

        #   solve
        self.analysis.solve()

        #   add to queue
        self.queue.put(self.analysis._solution_data)

    def _set_variables(self):
        """
        Function sets all values of defined variables of the MBD system
        :return:
        """
        for variable in self.MBD_system.variables:
            variable.set_value()
