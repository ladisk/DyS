'''
Created on 13. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''

import A_theta
import numpy as np


def dr_ij_P_dq(body_id_i, theta_i, u_iP, body_id_j, theta_j, u_jP):
    """
    Function calculates distance between two points on two different bodies in plane
    Args:
        R_i - vector of coordinates of centre of mass of body_i
        R_j - 
        theta_i - 
        theta_j -  
    Returns:
    matrix
    """
    if body_id_j == "ground":
        matrix_ = np.hstack((np.eye(2), np.array([np.dot(A_theta.A_theta_matrix(theta_i), u_iP)]).T, np.zeros([2, 3])))
    else:
        matrix_ = np.hstack((np.eye(2), np.array([np.dot(A_theta.A_theta_matrix(theta_i), u_iP)]).T, -np.eye(2), -np.array([np.dot(A_theta.A_theta_matrix(theta_j), u_jP)]).T))
    return matrix_