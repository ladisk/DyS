'''
Created on 27. apr. 2015

@author: lskrinjar

'''
import numpy as np


def get_dq_from_q(q):
    """
    Args:
        q - vector of position and velocity of all bodies of MBD system
        q = [q, dq]
    Returns:
        dq - vector of velocities - second half of vector q
    """
    dq = q[0.5 * len(q):]
    return dq


if __name__ == "__main__":
    _q = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    print get_dq_from_q(_q) 
    