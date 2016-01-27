"""
Created on 17. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""

import numpy as np


from A import A_matrix
from A import A_3D_matrix


def Ai_ui_P_vector(u_P, theta):
    """
    Function calculates vector of point P in LCS
    :param u_P:     a position vector or matrix of vectors (2D or 3D vectors)
    :param theta:   angle of rotation (in radians)
    :return vector: a position vector transformed with transformation matrix (rotation) - size is equal as u_P
                    position of point P in LCS of a body after rotation
    """
    #   check if array(s) is/are 2D(x, y) or 3D(x, y, z)
    u_P_dim = u_P.ndim

    #   if only one vector is passed as input parameter
    if u_P_dim == 1:
        n = len(u_P)
        #   2D(x, y)
        if n == 2:
            matrix = A_matrix(theta)

        #   3D(x, y, z)
        elif n == 3:
            matrix = A_3D_matrix(theta)

        else:
            pass

        vector = np.dot(matrix, u_P)
        return vector

    #   if an array of vector points (matrix) is passed as input parameter
    elif u_P_dim == 2:
        [rows, cols] = np.shape(u_P)
        #   predefine empty matrix
        matrix = np.zeros((rows, cols))

        #   2D(x, y)
        if cols == 2:
            _A_matrix = A_matrix(theta)

        #   3D(x, y, z)
        elif cols == 3:
            _A_matrix = A_3D_matrix(theta)
        else:
            pass

        for i, _u_P in enumerate(u_P):
            matrix[i, :] = np.dot(_A_matrix, _u_P)
        return matrix

    else:
        raise ValueError, "Wrong input!"


if __name__ == "__main__":

    #    vector 2D(x, y)
    # uP = np.random.rand(2)
    #   vector 3D(x, y, z)
    # uP = np.random.rand(3)
    #    matrix of 2D vectors
    # uP = np.random.rand(10, 2)
    #    matrix of 3D vectors
    uP = np.random.rand(10, 3)

    theta = 0
    print "input ="
    print uP
    _uP = Ai_ui_P_vector(uP, theta)
    print "output ="
    print _uP
