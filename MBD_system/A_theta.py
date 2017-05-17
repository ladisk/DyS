'''
Created on 7. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''

import numpy as np


def A_theta_matrix(theta):
    """
    First partial derivative of matrix to theta
    in:
        theta in radians
    out:
        transformation matrix
    """
    matrix_ = np.array([[-np.sin(theta), -np.cos(theta)],
                        [np.cos(theta), -np.sin(theta)]])
    
    return matrix_