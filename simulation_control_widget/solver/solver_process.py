"""

created by: lskrinjar
date of creation: 16/06/2016
time of creation: 17:57
"""
import os
import sys
import ctypes
import copy
import time
import multiprocessing
import traceback
import itertools
import threading
import itertools
import numpy as np

from kinematic_analysis import KinematicAnalysis
from integrate_euler import IntegrateEuler
from integrate_RKF45 import IntegrateRKF45


class SolverProcess(object):
    """

    """
    __id = itertools.count(0)

    def __init__(self, MBD_system):#MBD_system=None, queue=None, lock=None
        """

        :param MBD_system:
        :param queue:
        :param lock:
        """
        super(SolverProcess, self).__init__()

        #   error handling
        # self._pconn, self._cconn = multiprocessing.Pipe()
        #self.child_pipe, self.parent_pipe = multiprocessing.Pipe(duplex = True)
        self._exception = None

        #   name
        #self._name = multiprocessing.current_process().name

        #   solver process id
        #self.solver_process_id = self.__id.next()

        #   create folder to store computation data
        self.folder_path = "job_"+str(1).zfill(4)
        if not os.path.exists(self.folder_path):
            os.makedirs(self.folder_path)
        else:
            for root, dirs, files in os.walk(os.path.abspath(self.folder_path)):
                for f in files:
                    os.unlink(os.path.join(root, f))

        #   MBD system
        self.MBD_system = MBD_system

        #   queue
        self.queue = None

        #   lock
        self.lock = None

        #   deamon
        self.daemon = False

        #   solver (analysis)
        self.analysis = self._analysis_type()

        #   open GL properties
        self._update_display_type = "dt"
        self._update_display_types = ["step", "dt"]

        if self.analysis is not None:
            self._dt = self.MBD_system._dt

    #@property
    #def exception(self):
    #    if self.child_pipe.poll():
    #        self._exception = self.child_pipe.recv()
    #    return self._exception

    def _analysis_type(self):
        """

        :return:
        """
        analysis = None
        #   kinematic analysis
        if hasattr(self.MBD_system, "analysis_type"):
            if self.MBD_system.analysis_type == "kinematic":
                analysis = KinematicAnalysis(self.MBD_system, parent=self)

            #   dynamic analysis
            elif self.MBD_system.analysis_type == "dynamic":

                if self.MBD_system.integrationMethod == "euler":
                    analysis = IntegrateEuler(self.MBD_system, parent=self)

                elif self.MBD_system.integrationMethod == "RKF45":
                    analysis = IntegrateRKF45(self.MBD_system, parent=self)
        else:
            print Warning, "Analysis type not defined!"

        return analysis

    def run(self, i):
        """

        :return:
        """
        print "pid =", os.getpid(), "i =", i, "name =", multiprocessing.current_process().name
        try:
            for i in range(0, 100):
                print "pid =", os.getpid(), "i =", i
                sys.stdout.flush()
            self._cconn.send(None)

        except Exception as e:
            tb = traceback.format_exc()
            self._cconn.send((e, tb))

        # return i*i
        # return x*x
        # ctypes.windll.user32.MessageBoxA(0, str(i), "Sem v funkciji", 1)
        # print "i =", i
        # print "I'm running!"
        # with self.lock:
        #   update all variables
        # self._set_variables()

        #   solve
        # self.analysis.solve()
        # for i in range(0, 100000000):
        #     np.sin(i)
        #     # time.sleep(2)
        #
        #     print "id =", id(self), i
        #     sys.stdout.flush()
            # time.sleep(1)
        #   add to queue
        # self.queue.put(self.analysis._solution_data)

    def _set_variables(self):
        """
        Function sets all values of defined variables of the MBD system
        :return:
        """
        for variable in self.MBD_system.variables:
            variable.set_value()
