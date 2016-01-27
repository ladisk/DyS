'''
Created on 16. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''

import numpy as np


try:
    from A_theta import A_theta_matrix
except:
    None

def I_A_uP_i_matrix(u_P, theta):
    """
    Function construct a matrix = [I A_theta_i*u_Pi]
    Args:
        u_P - a position vector of point P in CM LCS
        theta in radians- angle of rotation of LCS
    Returns:
        matrix of shape (6, 2)
    """
    matrix_ =  np.hstack((np.eye(2), np.array([np.dot(A_theta_matrix(theta), u_P)]).T))
    return matrix_

if __name__ == "__main__":
    a = I_A_uP_i_matrix(np.random.rand(2), 1)
    print a
    