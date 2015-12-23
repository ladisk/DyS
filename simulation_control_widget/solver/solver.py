'''
Created on 10. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''

from pprint import pprint
import sys
import time

from PyQt4 import QtCore
from PyQt4.QtCore import *

import numpy as np
from solve_ODE import SolveODE


class FinishedSignal(QtCore.QObject):
    signal_finished = QtCore.pyqtSignal()
    
class runningSignal(QtCore.QObject):
    signal_running = QtCore.pyqtSignal(["QString"])

class stoppedSignal(QtCore.QObject):
    signal_stopped = QtCore.pyqtSignal(["QString"])

class Solver(QtCore.QThread):
    """
    Solver thread as object
    """
    def __init__(self, MBD_system, parent=None):
        QThread.__init__(self)
        self._parent = parent
        self.running = False
        self.running_signal = runningSignal()
        self.stopped = False
        self.stopped_signal = stoppedSignal()
        self.finished_signal = FinishedSignal()
        self.paused = False
         
        self.MBD_system = MBD_system
        
        self.solveODE = SolveODE(MBD_system=self.MBD_system, parent=self._parent)

    def start_solver(self):
        """
        Method that starts the ODE solver (integrtor)
        """
        self.running_signal.signal_running.emit("Running")

        while self.solveODE.running and not self.solveODE.stopped:
            self.solveODE.solve_ODE()

    def stop_solver(self):
        self.stopped_signal.signal_stopped.emit("Stopped")
