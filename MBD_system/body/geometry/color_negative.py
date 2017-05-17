"""

created by: lskrinjar
date of creation: 13/04/2016
time of creation: 21:11
"""
import numpy as np

def color_negative(color):
    """

    :param color:
    :return:
    """
    color = np.ones(3, dtype="float32") - color

    return color

if __name__ == "__main__":
    c1 = np.random.rand(3)
    print c1
    print color_negative(c1)