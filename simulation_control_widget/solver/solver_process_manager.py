"""
Created on 16. jun. 2016

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
# coding=utf-8
import sys
import os
import copy
import Queue
import threading
import multiprocessing
import traceback
from multiprocessing import Manager, Pool, Process, Lock, Queue
from multiprocessing.managers import BaseManager
from multiprocessing.pool import ThreadPool
from PyQt4 import QtCore


from simulation_control_widget.solver.solver import Solver
from simulation_control_widget.solver.solver_thread import SolverThread
from simulation_control_widget.solver.solver_process import SolverProcess


def f(mbd_system):
    print 1
    #caller = sys._getframe(2)
    #__locals = caller.f_locals
    #print 2, __locals
    #   evaluated force
    #my_MBD_system = __locals["MBD_system"]
    my_MBD_system = mbd_system
    print 3
    print my_MBD_system
    print 4
    sp = SolverProcess(my_MBD_system)
    print 5
    sp.run(2)
    print 6
    print x
    print 7
    # for i in range(0, 100):
    #     print "pid =", os.getpid(), "i =", i
    sys.stdout.flush()
    # return obj


class SolverProcessManager(QtCore.QThread):
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
        self.MBD_system.dys = None
        #self.MBD_system_0 = copy.copy(MBD_system)
        #self.MBD_system_1 = copy.copy(MBD_system)
        #self.MBD_system_2 = copy.copy(MBD_system)
        #self.MBD_system_3 = copy.copy(MBD_system)

        #   queue
        self.queue = Queue()

        #   lock
        #self.lock = threading.Lock()

        #   number of cpu
        self.cpu_num = multiprocessing.cpu_count()
        print "self.cpu_num =", self.cpu_num

        #   pool
        self.pool = Pool(processes=self.cpu_num-1)              # start 1 process per CPU core

        self.N = int(self._parent.MBD_system.monte_carlo_number_of_simulations)
        #   status
        self.running = False

        #   create processes
        self.processes = []

    def start_solver(self):
        """

        :return:
        """
        #   create list of MBD object with different values of selected variables
        self.N = 10
        # pool = Pool(processes=3)
        #for i in range(0, self.N):#self.N

        #    solver_process = SolverProcess(self.MBD_system)
            # a = pool.apply_async(f(i, solver_process), args=(i,solver_process,))
            # a.get(timeout=1)
        #    self.processes.append(solver_process)
        # pool.close()
        # pool.join()
        # print self.processes
        # self.processes[0].analysis._active_scene_GL = True
        self.running = True
        # [self.processes[i].runn(i) for i in range(0, 3)]
        results = []
        # for i in range(0, 10):
        #     a = self.pool.apply_async(self.processes[i].runn(i), args=(i,))
        #     results.append(a)
        #for i in range(0, self.N):
        input_data = [copy.copy(self.MBD_system) for i in range(3)]
        #input_data = [3 for i in range(3)]
        #MBD_system = self.MBD_system
        #self.pool.apply_async(f, args=(i, MBD_system, ))
        self.pool.map(f, input_data)

        self.pool.close()
        self.pool.join()

        for p in self.processes:
            if p.exception:
                error, traceback = p.exception
                print traceback

        # pool = Pool(processes=4)
        # a0 = pool.apply_async(self.processes[0].runn(i), args=(0,))
        # a1 = pool.apply_async(self.processes[1].runn(i), args=(1,))
        # a2 = pool.apply_async(self.processes[2].runn(i), args=(2,))
        # a3 = pool.apply_async(self.processes[3].runn(i), args=(3,))

        # for i in range(0, self.N):
        #     print i
        #     self.pool.apply_async(self.processes[i].runn(i), args=(i,))
            # self.processes[i].runn(i)
            # pid = os.spawnlp(os.P_NOWAIT, "notepad.exe", "")
            # self.pool.apply_async(self.processes[i].runn, args=(i,))
        # self.pool.map(1)
        # for i, process in enumerate(self.processes):
        #     self.pool.apply_async(self.processes[i])

        # self.pool.close()
        # self.pool.join()
            # print "process id =", id(process)
            #   start the process
            # process.start()

        # for process in self.processes:
        #     process.join()

        # for process in self.processes:
        #     process.join()
        #     print self.queue.get()

    def stop_solver(self):
        """

        :return:
        """
