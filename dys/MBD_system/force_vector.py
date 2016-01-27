'''
Created on 4. jun. 2014

@author: luka.skrinjar
'''
import itertools

import numpy as np


class Force_Q_e_vector(object):
    '''
    classdocs
    '''
    __id = itertools.cycle([0, 1])

    def __init__(self, f_s = [], force_Q_e_body_matrix = [], I_r_ij = []):
        '''
        Constructor
        Args:
        f_s - spring force - scalar value
        force_Q_e_body_matrix - 
        I_r_ij - unit vector in direction of a spring
        
        Returns:
        
        '''
        self.id = self.__id.next()
        
        if f_s == [] and force_Q_e_body_matrix == [] and I_r_ij == []:
            self.vector = np.zeros(3)
        else:
            
            if self.id == 0:
                # print "I_r_ij =", I_r_ij
                # print "force_Q_e_body_matrix =", force_Q_e_body_matrix
                # print "np.dot(force_Q_e_body_matrix, I_r_ij) =", np.dot(force_Q_e_body_matrix, I_r_ij)
                self.vector = -f_s * np.dot(force_Q_e_body_matrix, I_r_ij)
                # print "self.vector =", self.vector
            elif self.id == 1:
                self.vector = +f_s * np.dot(force_Q_e_body_matrix, I_r_ij)
            else:
                raise ValueError, "Wrong id!"                                  