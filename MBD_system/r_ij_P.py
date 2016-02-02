"""
Created on 13. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""

import A
import numpy as np


def r_ij_P(R_i, theta_i, u_iP, R_j, theta_j, u_jP):
    """
    Function calculates distance between two points on two different bodies in plane
    Args:
        R_i - vector of coordinates of centre of mass of body_i
        R_j - 
        theta_i - 
        theta_j -  
    """
    vector = R_i + np.dot(A.A_matrix(theta_i), u_iP) - R_j - np.dot(A.A_matrix(theta_j), u_jP)
    return vector