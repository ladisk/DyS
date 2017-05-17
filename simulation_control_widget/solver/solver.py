"""
Created on 10. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import os
from multiprocessing import Manager, Pool, Process
from PyQt4 import QtCore


from dynamic_analysis import DynamicAnalysis
from kinematic_analysis import KinematicAnalysis
from integrate_euler import IntegrateEuler
from integrate_RKF45 import IntegrateRKF45
from integrate_HHT_I3 import IntegrateHHT_I3


class Solver(QtCore.QThread):#QRunnable, QRunnable
    """
    Solver thread as object
    """
    def __init__(self, MBD_system, parent=None):
        QtCore.QThread.__init__(self)
        # super(Solver, self).__init__(parent)
        self._parent = parent
        self.running = False
        self.stopped = False
        self.paused = False

        #   pointer to MBD system
        self.MBD_system = MBD_system

        #   solver (analysis)
        self.analysis = self._analysis_type()
        self.analysis._active_scene_GL = True

        #   process id
        self.pid = os.getpid()

    def _analysis_type(self):
        """
        Select analysis type
        :return:
        """
        #   kinematic analysis
        if self.MBD_system.analysis_type == "kinematic":
            analysis = KinematicAnalysis(self.MBD_system, parent=self._parent)

        #   dynamic analysis
        elif self.MBD_system.analysis_type == "dynamic":
            if self.MBD_system.integrationMethod.lower() == "euler":
                analysis = IntegrateEuler(self.MBD_system, parent=self._parent)
            
            elif self.MBD_system.integrationMethod == "RKF45":
                analysis = IntegrateRKF45(self.MBD_system, parent=self._parent)

            elif self.MBD_system.integrationMethod == "HHT-I3":
                analysis = IntegrateHHT_I3(self.MBD_system, parent=self._parent)

        else:
            analysis = None
            raise ValueError, "Analysis type not defined!"

        return analysis

    def start_solver(self):
        """
        Method that starts the ODE solver (integrtor)
        """
        #   solver
        # self.analysis = self._analysis_type()
        if self.analysis is not None:
            self.analysis.solve()
        else:
            raise ValueError, "Analysis type not specified by user!"

    def stop_solver(self):
        """

        :return:
        """
