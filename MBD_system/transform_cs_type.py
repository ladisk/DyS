# coding=utf-8

import numpy as np


def transform_cartesian2polar(x, y):
    """

    :param x:
    :param y:
    :return:
    """
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi


def transform_polar2cartesian(r, phi):
    """

    :param r:
    :param phi:
    :return:
    """
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return np.array([x, y])


if __name__ == "__main__":
    x = 1.
    y = 0.

    print transform_cartesian2polar(x, y)