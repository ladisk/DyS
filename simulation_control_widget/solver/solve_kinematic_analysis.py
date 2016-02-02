"""
Created on 29. jan. 2016

@author: luka.skrinjar
"""
from pprint import pprint
import scipy.optimize
from PyQt4.QtGui import *
from PyQt4.QtCore import *


from DAE_fun import DAEfun
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2R_i import q2R_i
from MBD_system import convert_bytes_to_
from MBD_system.check_filename import check_filename
from MBD_system.solution_data.solution_data import SolutionData
from solve_dynamic_analysis import SolveDynamicAnalysis


class SolveKinematicAnalysis(SolveDynamicAnalysis):
    """
    classdocs
    """
    def __init__(self, MBD_system, parent=None):
        """
        Constructor
        """
        super(SolveKinematicAnalysis, self).__init__(MBD_system, parent)
        #   parent
        self._parent = parent

        #    DAE fun object
        self.DAE_fun = DAEfun(self._MBD_system, parent=self)

    def solve(self):
        """
        Solves a kinematic system
        """
        print "solve()"
        # pprint(vars(self))

        #   start solve
        self.start_solve()


        self.simulation_id = 0
        self.FLAG = 1


        # scipy.optimize.newton()

        t = 0
        h = self._MBD_system.Hmax
        self.t_n = self._MBD_system.t_n

        #   initial approximation
        q = self.DAE_fun.evaluate_q0()

        #   kinematic analysis
        while self.FLAG == 1:
            print "t =", t
            if self.stopped:
                # self.update_GL_(t=t, q=w)
                self.stop_solve()
                self.FLAG = 0

            if t >= self.t_n:
                self.FLAG = 0
                self.finished = True
                print "finished!"
                # self.update_GL_(t=t, q=w)

            if self.finished or self.failed:
                self.finished_solve()


            #   evaluate C vector
            C = self.DAE_fun.evaluate_C(t, q)
            print "C =", C


            #   increase time step
            self.t = t = t + h