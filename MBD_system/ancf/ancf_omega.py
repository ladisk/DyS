"""

created by: lskrinjar
date of creation: 20/03/2017
time of creation: 22:10
"""
import numpy as np


def ancf_omega(e_x, de_x):
    """

    :return:
    """
    omega = (e_x[1] * de_x[0] - e_x[0] * de_x[1])/(np.linalg.norm(e_x) ** 2)

    return omega