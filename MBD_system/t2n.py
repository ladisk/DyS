"""

created by: lskrinjar
date of creation: 23/01/2016
time of creation: 15:00
"""

import numpy as np

def t2n(t):
    """
    Function evaluates normal based on tanget vector
    :param t:   tangent
    :return n:  normal
    """
    n = np.array([t[1], -t[0]])

    return n
