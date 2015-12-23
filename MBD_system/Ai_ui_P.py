'''
Created on 17. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''

import numpy as np


from A import A_matrix
from A import A_3D_matrix


def Ai_ui_P_vector(u_P, theta):
    """
    Function calculates vector of point P in LCS
    Args:
        u_P (in CM LCS) - 
        theta - angle of rotation (in radians)
    Returns:
        vector_ - position of point P in CS LCS of a body
    """
    if len(u_P) == 2:
        matrix_ = A_matrix(theta)
    if len(u_P) == 3:
        matrix_ = A_3D_matrix(theta)
    #    if u_P is only a vector of one point
    vector_ = np.dot(matrix_, u_P)
    return vector_
    #
    # #    if u_P is a matrix of points
    # else:
    #     matrix_ = np.zeros(np.shape(u_P))
    #     for i, u_P_i in enumerate(u_P):
    #         matrix_[i, :] = np.dot(A_matrix_, u_P_i)
    #     return matrix_


if __name__ == "__main__":
    #    2D array
    # print Ai_ui_P_vector(np.random.rand(2), 1)
    #    3D array - NOT WORKING
    print Ai_ui_P_vector(np.random.rand(3), 1)
    #    2D matrix
    # m = np.random.rand(10, 2)
    # print m

    