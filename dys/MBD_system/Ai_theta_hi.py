'''
Created on 17. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''

import numpy as np


try:
    from A_theta import A_theta_matrix
except:
    None

def Ai_theta_hi(h, theta):
    """
    Function calculates vector h_i
    Args:
        u_P (in CM LCS) - 
        theta - angle of rotation
    Returns:
        vector_ - position of point P in CS LCS of a body
    """

    vector_ = np.dot(A_theta_matrix(theta), h)

    return vector_

# if __name__ == "__main__":
#     print Ai_theta_hi(np.random.rand(2), 1)
#     