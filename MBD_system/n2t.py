"""

created by: lskrinjar
date of creation: 23/01/2016
time of creation: 14:52
"""
import numpy as np

def n2t(n):
    """
    Function evaluates tangent based on normal vector
    :param n:   normal
    :return t:  tangent
    """
    t = np.array([-n[1], n[0]])

    return t
