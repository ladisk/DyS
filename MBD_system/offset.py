"""

created by: lskrinjar
date of creation: 17/06/2016
time of creation: 18:24
"""
import numpy as np
from find_nearest import find_nearest

def offset(x, y, x0=None, xn=None, dy=None):
    """

    :param x:
    :param y:
    :param x0:
    :param xn:
    :return:
    """
    if xn is not None and xn is not None:
        dy = np.mean(y[find_nearest(x, x0):find_nearest(x, xn)])

    return y - dy