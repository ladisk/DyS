"""
Created on 16. jun. 2016

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
# coding=utf-8
import copy
import Queue
import threading
import multiprocessing
from multiprocessing import Manager, Pool, Process, Lock, Queue
from multiprocessing.managers import BaseManager
from multiprocessing.pool import ThreadPool
from PyQt4 import QtCore


from simulation_control_widget.solver.solver import Solver
from simulation_control_widget.solver.solver_thread import SolverThread

class SolverThreadManager(QtCore.QThread):
    """
    Manager of simulation threads
    """

    def __init__(self, MBD_system, parent=None):
        QtCore.QThread.__init__(self)
        #   parent
        self._parent = parent

        #   analysis
        self.analysis = None

        #   MBD system
        self.MBD_system = copy.copy(MBD_system)

        #   number of cpu
        self.cpu_num = multiprocessing.cpu_count()

        #   queue
        self.queue = Queue()

        #   lock
        self.lock = threading.Lock()

        #   create processes
        self.processes = []
        for i in range(0, int(self._parent.MBD_system.monte_carlo_number_of_simulations)):

            solver_process = SolverThread(MBD_system, self.queue, self.lock)

            self.processes.append(solver_process)

    def start_solver(self):
        """

        :return:
        """
        self.processes[0].analysis._active_scene_GL = True

        for process in self.processes:
            print "process id =", id(process)
            #   start the process
            process.start()

        for process in self.processes:
            process.join()

        # for process in self.processes:
        #     process.join()
        #     print self.queue.get()

    def stop_solver(self):
        """

        :return:
        """
