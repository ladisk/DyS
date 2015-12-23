'''
Created on 7. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''

import numpy as np


def A_matrix(theta):
    """

    :param theta: angle in rad
    :return:
    :rtype:
    >>> 1+1
    2
    >>>
    """
    matrix = np.array([[np.cos(theta), -np.sin(theta)],
                       [np.sin(theta), np.cos(theta)]])
    return matrix

def A_3D_matrix(theta):
    """

    :param theta: angle in rad in 3D with only one rotation about z axis
    :return:
    """
    matrix = np.array([[np.cos(theta), -np.sin(theta), 0],
                       [np.sin(theta), np.cos(theta), 0],
                       [0, 0, 1]])
    return matrix


if __name__ == "__main__":
    pass