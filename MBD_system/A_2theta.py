'''
Created on 2. maj 2015

@author: lskrinjar

'''

import numpy as np


def A_2theta_matrix(theta):
    """
    Second partial derivative of matrix to theta
    in:
        theta in radians
    out:
        transformation matrix
    """
    matrix_ = np.array([[-np.cos(theta), np.sin(theta)],
                        [-np.sin(theta), -np.cos(theta)]])
    return matrix_