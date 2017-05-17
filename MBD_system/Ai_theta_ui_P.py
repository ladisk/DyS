"""
Created on 17. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""

import numpy as np


try:
    from A_theta import A_theta_matrix
except:
    None

def Ai_theta_ui_P_vector(u_P, theta, _id):
    """
    Function calculates vector of point P in LCS based on formula A_theta(theta_i).ui_P
    Args:
        u_P (in CM LCS) - 
        theta - angle of rotation in radians
        id - i=0 or j=1
    Returns:
        vector_ - position of point P in CS LCS of a body
    """
    if _id == 0:
        vector = +np.dot(A_theta_matrix(theta), u_P)

    elif _id == 1:
        vector = -np.dot(A_theta_matrix(theta), u_P)

    else:
        raise ValueError, "id not correct!"

    return vector


if __name__ == "__main__":
    print Ai_theta_ui_P_vector(np.random.rand(2), 1)

