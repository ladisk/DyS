'''
Created on 7. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''

import numpy as np


def A_inv_matrix(theta):
    """
    
    """
    matrix = np.array([[np.cos(theta), np.sin(theta)],
                       [-np.sin(theta), np.cos(theta)]])
    
    return matrix

if __name__ == "__main__":
    pass