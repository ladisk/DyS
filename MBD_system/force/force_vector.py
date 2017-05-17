"""
Created on 4. jun. 2014

@author: luka.skrinjar
"""
import itertools


import numpy as np


from global_variables import GlobalVariables


class Force_Q_e_vector(object):
    """
    classdocs
    """
    __id = itertools.cycle([0, 1])

    def __init__(self, f_s=0., force_Q_e_body_matrix=[], I_r_ij=[], n=3):
        """
        Constructor
        Args:
        f_s - spring force - scalar value
        force_Q_e_body_matrix - 
        I_r_ij - unit vector in direction of a spring
        n - size of force vector (planar)
            for rigid body is 3
            for flexible body is equal to number of absolute nodal coordinated of a mesh
        """

        self.id = self.__id.next()

        #   size of force vector
        self.n = n
        
        if f_s == 0. and force_Q_e_body_matrix == [] and I_r_ij == []:
            self.Q_e = np.zeros(self.n)
        else:
            
            if self.id == 0:
                self.Q_e = -f_s * np.dot(force_Q_e_body_matrix, I_r_ij)

            elif self.id == 1:
                self.Q_e = +f_s * np.dot(force_Q_e_body_matrix, I_r_ij)

            else:
                raise ValueError, "Wrong id!"

    def evaluate_Q_e(self, M=0.):
        """

        :return:
        """
        self.Q_e[-1] = M
