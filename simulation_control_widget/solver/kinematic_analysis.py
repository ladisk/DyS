"""
Created on 29. jan. 2016

@author: luka.skrinjar
"""
import numpy as np
import scipy

from DAE_fun import DAEfun
from dynamic_analysis import DynamicAnalysis


class KinematicAnalysis(DynamicAnalysis):
    """
    classdocs
    """
    def __init__(self, MBD_system, parent=None):
        """
        Constructor
        """
        super(KinematicAnalysis, self).__init__(MBD_system, parent)
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

        t = 0.
        h = self._MBD_system.Hmax
        self.t_n = self.MBD_system.t_n

        #   initial approximation
        q = self.q_sol = self.MBD_system.evaluate_q0()
        #   use  positions vector (translational and rotational)
        q = q[0:len(q)/2]
        #   us
        #   solution containers
        self._solution_containers()

        #   size of C, C_q, C_t vectors, matrix
        self.DAE_fun._C_q_size()

        #   options to track data
        self.q_sol_matrix = np.array([q])
        self._energy_t = 0
        self._energy_delta = np.nan

        #   check before run
        if not self.DAE_fun.check_C(q, t):
            self.FLAG = 1
        else:
            print "Error found before simulation run!"
            self.FLAG = 0
            self.finished_solve()

        #   kinematic analysis
        while self.FLAG == 1:
            # print "----------------------------------------"
            # print "t =", t
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
            # C = self.DAE_fun.evaluate_C(q, t)
            # print "C =", len(C)
            # print C


            # C_q = self.DAE_fun.evaluate_C_q(q)
            # print "C_q ="
            # print C_q
            # print "q (in) =", q
            q = scipy.optimize.fsolve(self.DAE_fun.evaluate_C, q, args=(t), fprime=self.DAE_fun.evaluate_C_q, xtol=self.MBD_system.TOL_dq_i, maxfev=50)
            # q = scipy.optimize.newton(self.DAE_fun.evaluate_C, q, fprime=self.DAE_fun.evaluate_C_q, tol=self.MBD_system.TOL_dq_i, maxiter=50)
            # print "q (out) =", q

            #   evaluate C_q matrix
            C_q = self.DAE_fun.evaluate_C_q(q, t)
            # print "C_q ="
            # print C_q
            #   evaluate C_t vector
            C_t = self.DAE_fun.evaluate_C_t(q, t)
            # print "C_t =", C_t
            #   solve linear system of equations for dq
            dq = np.linalg.solve(C_q, -C_t)

            #   assemble new q vector q=[q, dq]
            _q = np.concatenate([q, dq])
            self.q_sol = _q
            #   evaluate vector Q_d
            Q_d = self.DAE_fun.evaluate_Q_d(_q)

            #   solve linear system of equations for ddq
            ddq = np.linalg.solve(C_q, Q_d)

            #   track data
            self.R = np.nan
            self._track_data(h, t, q)

            self._info(t, q)

            #   increase time step
            t = t + h
            self.t = t
        
        #    integration finished - end time reached successfully
        if self.finished or self.failed:
            self.finished_solve()
        
        #    save data to file
        self.write_simulation_solution_to_file()