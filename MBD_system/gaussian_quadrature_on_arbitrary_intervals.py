"""

created by: lskrinjar
date of creation: 09/08/2016
time of creation: 20:16
"""

def change_interval(x, a, b):
    """
    Function changes value of x that is is inside interval [-1, +1]
    According to Numerical Analysis 9th ed - Burden and Faires, 2010
    :param a:
    :param b:
    :return:
    """
    x1 = ((b - a) * x + (b + a)) / 2.

    interval_ab = ((b - a) / 2.)

    return x1, interval_ab