'''
Created on 30. maj 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''
from A import A_matrix
import numpy as np


def hi_Ai_theta_ui_P_constant(h, theta, u_P):
    """
    Function calculates one component for C_q matrix for revolute joint
    Args:
    h - 
    theta - 
    u_P - 
    """
    constant = np.dot(h, np.dot(A_matrix(theta), u_P))
    return constant


if __name__ == "__main__":
    h = np.array([1, 2])
    theta = 0.1
    u_P = np.array([1, 10])
    print hi_Ai_theta_ui_P_constant(h, theta, u_P)