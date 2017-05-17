'''
Created on 30. maj 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''
from A_theta import A_theta_matrix
import numpy as np
from r_ij_P import r_ij_P


def r_ij_P_Ai_theta_hi_constant(rijP, theta_i, hi):
    """
    
    Args:
    rijP
    hi
    """
    
    Ai_ = A_theta_matrix(theta_i)
    
    constant = np.dot(rijP, np.dot(Ai_, hi))
    
    return constant


if __name__ == "__main__":
    R_i = np.array([1, 1])
    theta_i = 0.1
    u_iP = np.array([1, 1])
    R_j = np.array([1, 1])
    theta_j = 0.2
    u_jP = np.array([-1, 2])
    hi = np.array([1, 10])
    print rijP_Ai_theta_hi_constant(R_i, theta_i, u_iP, R_j, theta_j, u_jP, hi)
