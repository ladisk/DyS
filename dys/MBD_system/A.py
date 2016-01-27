'''
Created on 7. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''

import numpy as np


def A_matrix(theta):
    """

    :param theta: angle in rad
    :return:
    """
    matrix = np.array([[np.cos(theta), -np.sin(theta)],
                       [np.sin(theta), np.cos(theta)]])
    
    return matrix

if __name__ == "__main__":
    pass