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

    def solve(self):
        """
        Solves a kinematic system
        """
        print "solve()"
        pprint(vars(self))

        #   start solve
        self.start_solve()
        # scipy.optimize.newton()